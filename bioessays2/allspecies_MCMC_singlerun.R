# clear the workspace
rm(list = ls(all = T))

library(rjags)
library(R2OpenBUGS)
library(coda)
library(extraDistr)
library(ggplot2)

##### specify Ben's post.summ function
post.summ = function(post, var) {
  
  # coerce to matrix for easy subsetting
  post.samp = as.matrix(post)
  
  # if parameter is indexed
  if(substr(var, nchar(var), nchar(var)) == "[") {
    # extract columns with headers equal to the desired variable
    post.sub = post.samp[,substr(colnames(post.samp), 1, nchar(var)) == var]
    # calculate desired quantities
    summ = apply(post.sub, 2, function(x) c(mean = mean(x), sd = sd(x), quantile(x, c(0.5, 0.025, 0.975))))
    return(summ)
  }
  
  # if parameter is not indexed
  if(substr(var, nchar(var), nchar(var)) != "[") {
    # extract the column with the same header as the desired variable
    post.sub = post.samp[,substr(colnames(post.samp), 1, nchar(var)) == var]
    # calculate the desired quantities
    summ = c(mean = mean(post.sub), sd = sd(post.sub), quantile(post.sub, c(0.5, 0.025, 0.975)))
    return(summ)
  }
}

#########################################################
#### ESTIMATION OF Q (RATE OF SUCCESS OF EACH SIRE) 
#### AND CALCULATION OF M (# MATES) FROM R (# SIRES)
#########################################################

##### read and prepare data #####
mamm <- read.csv("/Users/Hannah/Dropbox/Paternity2/BayesianMCMC/mammals_Avise/single_run/paternity_mammals_Avise.csv")
fish <- read.csv("/Users/Hannah/Dropbox/Paternity2/BayesianMCMC/fish/single_run/paternity_fish.csv")
herp <- read.csv("/Users/Hannah/Dropbox/Paternity2/BayesianMCMC/herps/single_run/paternity_herps.csv")
invert <- read.csv("/Users/Hannah/Dropbox/Paternity2/BayesianMCMC/inverts/single_run/paternity_inverts.csv")
dat <- rbind(mamm, fish, herp, invert)

# expand data to create multiple broods within each species (sample size * 10)
for(i in 1:nrow(dat)){
  temp.sp <- expand.grid(species = rep(dat[i,]$species, dat[i,]$nbrood*10))
  if(i == 1){
    temp <- temp.sp
  } else {
    temp <- rbind(temp, temp.sp)
  } 
}
mates.df <- merge(temp, dat, all.x = TRUE, all.y = FALSE)


##### specify model code #####
model_byspecies_q.file <- "model_allspecies_q.txt"
jagsscript.byspecies_q <- cat("
  model{
  # define prior: probability of q (prob of success for each sire)
  q ~ dbeta(1,1)
  
  # likelihood: stochastic relationship
  # sires ~ (truncated binomial)
  for(i in 1:N){
    nsires[i] ~ dbinom(q, littersize[i]) T(0,)
  }
  
  # DERIVED QUANTITIES
  # Calculate the number of mates (deterministic relation?)
  # q is probability of success here
  # for(i in 1:N){
  #   nmates[i] <- (q * nsires[i]) / (1-q)
  # }
  
  # mates ~ (negative binomial); 
  ##
  # NOTE: in R, representation of NB dist is same as Casella-Berger stat theory book
  # whereas, Ash uses the alternative NB pmf seen in the note for NB dist in Casella-Berger
  # Therefore, P(X=x|r,p) used in R is P(X=M-r|r,q) in terms of number of failures M-r
  # mean of NB is μ = n(1-p)/p and variance n(1-p)/p^2, where p is the prob of SUCCESS
  # R is calculating the number of failures M-r, and we generated the number of successes (sires, R)
  # therefore the number of mates M = failures + successes = failed.m + nsires
  ##
  for(i in 1:N){
    failed.m[i] ~ dnbinom(q, nsires[i]) T(nsires[i],)
    nmates[i] <- failed.m[i] + nsires[i]
  }
  for(i in 1:N){
    est.M[i] <- nsires[i]/q
  }
  
  # Calculated quantities from original data
  for(j in 1:J){
    avg.R[j] <- avgbrood[j]*q/(1-(1-q)^avgbrood[j])     # avg sires, given avgbrood and estimated q
  }
  for(j in 1:J){
    avg.M[j] <- avgsire[j]/q     # avg mates, given avgsire and estimated q
  }
  for(j in 1:J){
    est.pmult.kq[j] <- (1 - (avgbrood[j] * q) * (1 - q)^(avgbrood[j] - 1)) / (1 - (1 - q)^(avgbrood[j]))   
    # estimated prob. of mult. paternity, given avgbrood and estimated q
  }
}
",file = model_byspecies_q.file)



#### BEGINNING OF REGENERATION OF DATA LOOP ####
# Remove any datapoints with NA for avgsire and avgbrood, since Poisson data gen will fail
mates.dat <- mates.df[!is.na(mates.df$avgsire),]

set.seed(36849)

#### (1) generate broods and sires from Poisson for each species ####
mates <- cbind(mates.dat, littersize = NA, nsires = NA)
for(i in 1:nrow(mates)){
  # generate broods
  nsample <- mates[i,]$nbrood   # number of broods
  avg.bs <- mates[i,]$avgbrood  # average brood size
  gen.k <- rtpois(1, avg.bs, a = 0, b = Inf)
  mates[i,]$littersize <- gen.k
  # generate sires
  avg.sire <- mates[i,]$avgsire  # average number of sires
  sire.max <- mates[i,]$littersize # max number of sires per litter is litter size
  gen.r <- rtpois(1, avg.sire, a = 0, b = sire.max)
  mates[i,]$nsires <- gen.r
}
## remove multiple paternity proportion of 1 or 0
#mates <- mates[mates$pmult != 0 & mates$pmult != 1,]
## remove an average sire value of NA
mates <- mates[!is.na(mates$avgsire), ]

# ##### (2) initial values for Q #####
# # Calculate an approximate initial q based on p (multiple paternity)
# ## Negative log-likelihood of zero-trunc binomial
# truncbin.lik <- function(q, s, k){
#   logl <- sum(s)*log(q) + sum(k-s)*(log(1-q)) - sum(log(1 - (1-q)^k))
#   return(-logl)
# }
# fq <- function(q, k, p) 1 - k*q*(1-q)^(k-1)/(1-(1-q)^k) - p
# mates.p <- mates
# mates.p$q0 <- mates.p$qs <- NA
# for(i in sort(unique(mates.p$species))){
#   for(j in 1:nrow(mates.p)){
#     mates.p[j,]$q0 <- uniroot(fq, c(0.01,.99), k=mates.p[j,]$littersize, p=mates.p[j,]$pmult)$root    # estimated prob success per sire
#   }
#   mates.p[mates.p$species==i,]$qs <- optim(.5, truncbin.lik, s=mates[mates.p$species==i,]$nsires, k=mates[mates.p$species==i,]$littersize, method = "Brent", hessian = TRUE, lower=0, upper=1)$par
# }
# q.init1 <- mean(mates.p$qs)
#   
# # obtain different q.init for 2nd chain
# q.init2 <- mean(sample(mates.p$qs, size = 25, replace = TRUE))
# #q.init2 <- mean(sample(mates.p$qs, size = 10, replace = TRUE))
# #q.init2 <- mean(sample(mates.p$qs, size = 10, replace = FALSE))
# #q.init2 <- 0.5084281
# #q.init2 <- 0.6295649
#   
# inits <- list(list(q=q.init1), list(q=q.init2))

##### (3) MCMC dimensions #####
ni = 10000
nb = 1000
nc = 2 # needs to match number of initial values proposed
nt = 2
n.iter = ni + nb

##### (4) parameters to monitor #####
params = c("q", "nsires", "failed.m", "nmates", "est.M", "avg.R", "avg.M", "est.pmult.kq")
# params = c("q", "littersize", "nsires", "failed.m", "nmates", "est.M", "avg.R", "avg.M", "est.pmult.kq")


##### CREATE LIST TO HOLD ALL RESULTS FROM JAGS
# # max number of iterations through all combinations of sires and litters
# maxlitter.for.minsire <- length(min(mates.dat$nsires):max(mates.dat$littersize))
# maxlitter.for.maxsire <- length(max(mates.dat$nsires):max(mates.dat$littersize))
# nlist <- sum(maxlitter.for.minsire:maxlitter.for.maxsire) 
# # slightly more slots in list than necessary, e.g. sires 5 & 6 only have max of 13 littersize
# # could change with random generation of litter sizes though
jags.results <- list() 

##### (5) run the model in JAGS #####
# Initial conditions for loop
# Likelihood needs to change for each nsire value, so q calculated for each # nsire
mates.p <- mates
species <- sort(unique(mates.p$species))
b <- 1
# Start JAGS
#starttime <- Sys.time()
for(r in species){
  # data containing only relevant values
  m.dat <- mates.p[mates.p$species==r,]
  orig.dat <- dat[dat$species==r,]
  jags.dat <- list(nsires = m.dat$nsires, 
                   littersize = m.dat$littersize, 
                   N = nrow(m.dat), 
                   avgbrood = orig.dat$avgbrood,
                   avgsire = orig.dat$avgsire,
                   J = nrow(orig.dat))
  
  # run JAGS (not using initial values for q)
  jmod <- jags.model(file = model_byspecies_q.file, data = jags.dat, n.chains = nc, n.adapt = 1000)  
  update(jmod, n.iter = nb, by = 1, progress.bar = 'text')
  post <- coda.samples(jmod, params, n.iter = ni, thin = nt) 
  
  # save current JAGS output
  name <- paste0("species=",r)
  jags.results[[name]] <- post
  
  # Save quantities from data
  MCMC_temp <- cbind(
    as.character(orig.dat$species),
    as.numeric(orig.dat$avgbrood),
    as.numeric(orig.dat$avgsire),
    as.numeric(orig.dat$pmult),
    post.summ(post, "q")[[1]],  # mean of the est. param. q
    post.summ(post, "q")[[2]],  # sd of the est. param. q
    post.summ(post, "q")[[4]],  # 2.5% est. param. q
    post.summ(post, "q")[[5]],  # 97.5% est. param. q
    post.summ(post, "nsires")[[1]],  # mean of nsires
    post.summ(post, "nsires")[[2]],  # sd of nsires
    post.summ(post, "nsires")[[4]],  # 2.5% nsires
    post.summ(post, "nsires")[[5]],  # 97.5% nsires
    post.summ(post, "nmates")[[1]],  # mean of nmates
    post.summ(post, "nmates")[[2]],  # sd of nmates
    post.summ(post, "nmates")[[4]],  # 2.5% nmates
    post.summ(post, "nmates")[[5]],  # 97.5% nmates
    post.summ(post, "est.M")[[1]],  # mean of est.M
    post.summ(post, "est.M")[[2]],  # sd of est.M
    post.summ(post, "est.M")[[4]],  # 2.5% est.M
    post.summ(post, "est.M")[[5]],  # 97.5% est.M
    post.summ(post, "est.pmult.kq")[[1]],  # mean of the est.pmult.kq
    post.summ(post, "est.pmult.kq")[[2]],  # sd of the est.pmult.kq
    post.summ(post, "est.pmult.kq")[[4]],  # 2.5% est.pmult.kq
    post.summ(post, "est.pmult.kq")[[5]],  # 97.5% est.pmult.kq
    post.summ(post, "avg.R")[[1]],  # mean of avg.R
    post.summ(post, "avg.R")[[2]],  # sd of avg.R
    post.summ(post, "avg.R")[[4]],  # 2.5% avg.R
    post.summ(post, "avg.R")[[5]],  # 97.5% avg.R
    post.summ(post, "avg.M")[[1]],  # mean of avg.M
    post.summ(post, "avg.M")[[2]],  # sd of avg.M
    post.summ(post, "avg.M")[[4]],  # 2.5% avg.M
    post.summ(post, "avg.M")[[5]]  # 97.5% avg.M
  )
  
  if(b==1) MCMC_sumtemp <- MCMC_temp else MCMC_sumtemp <- rbind(MCMC_sumtemp, MCMC_temp)
  
  # move to next step
  b <- b + 1
}

MCMC_summary <- data.frame(species = as.character(MCMC_sumtemp[,1]),
                           avgbrood = as.numeric(MCMC_sumtemp[,2]),
                           avgsire = as.numeric(MCMC_sumtemp[,3]),
                           pmult = as.numeric(MCMC_sumtemp[,4]),
                           mean.q = as.numeric(MCMC_sumtemp[,5]),
                           sd.q = as.numeric(MCMC_sumtemp[,6]),
                           q2.5 = as.numeric(MCMC_sumtemp[,7]),
                           q97.5 = as.numeric(MCMC_sumtemp[,8]),
                           mean.nsires = as.numeric(MCMC_sumtemp[,9]),
                           sd.nsires = as.numeric(MCMC_sumtemp[,10]),
                           nsires2.5 = as.numeric(MCMC_sumtemp[,11]),
                           nsires97.5 = as.numeric(MCMC_sumtemp[,12]),
                           mean.nmates = as.numeric(MCMC_sumtemp[,13]),
                           sd.nmates = as.numeric(MCMC_sumtemp[,14]),
                           nmates2.5 = as.numeric(MCMC_sumtemp[,15]),
                           nmates97.5 = as.numeric(MCMC_sumtemp[,16]),
                           mean.estM = as.numeric(MCMC_sumtemp[,17]),
                           sd.estM = as.numeric(MCMC_sumtemp[,18]),
                           estM2.5 = as.numeric(MCMC_sumtemp[,19]),
                           estM97.5 = as.numeric(MCMC_sumtemp[,20]),
                           mean.est.pmult = as.numeric(MCMC_sumtemp[,21]),
                           sd.est.pmult = as.numeric(MCMC_sumtemp[,22]),
                           est.pmult2.5 = as.numeric(MCMC_sumtemp[,23]),
                           est.pmult97.5 = as.numeric(MCMC_sumtemp[,24]),
                           mean.avgR = as.numeric(MCMC_sumtemp[,25]),
                           sd.avgR = as.numeric(MCMC_sumtemp[,26]),
                           avgR2.5 = as.numeric(MCMC_sumtemp[,27]),
                           avgR97.5 = as.numeric(MCMC_sumtemp[,28]),
                           mean.avgM = as.numeric(MCMC_sumtemp[,29]),
                           sd.avgM = as.numeric(MCMC_sumtemp[,30]),
                           avgM2.5 = as.numeric(MCMC_sumtemp[,31]),
                           avgM97.5 = as.numeric(MCMC_sumtemp[,32]))

save(MCMC_summary, file = paste0("MCMC_summary_singlerun.rda"))


# plot single run with original pmult data (from Avis)
pdf(file = paste0("p_singlerun_log.pdf"), width = 11, height = 7.5)
print(
  ggplot() +
    stat_smooth(data = MCMC_summary, aes(log(avgbrood), mean.est.pmult), color = "red", method = "loess") +
    geom_point(data = MCMC_summary, aes(log(avgbrood), mean.est.pmult), color = "red", alpha = 0.2, size = 2.5) +
    geom_point(data = dat, aes(log(avgbrood), pmult)) +
    #stat_smooth(data = dat, aes(log(avgbrood), pmult), color = "red") +
    labs(x="Log litter/clutch size", y="Probability of multiple paternity") +
    lims(y = c(0,1)) +
    scale_x_continuous(breaks = c(seq(0,15, by = 2)) , limits = c(0,15)) +
    theme_bw(base_size = 16)
)
dev.off()


## Color by taxa group
MCMC_sum2 <- MCMC_summary
MCMC_sum2$group <- c(rep("mammals", 49), rep("fish", 28), rep("herps", 16), rep("inverts", 30))
dat2 <- dat
dat2$group <- c(rep("mammals", 49), rep("fish", 28), rep("herps", 18), rep("inverts", 30))

pdf(file = paste0("p_singlerun_group.pdf"), width = 11, height = 7.5)
print(
  ggplot() +
    stat_smooth(data = MCMC_sum2, aes(log(avgbrood), mean.est.pmult), color = "red", method = "loess") +
    geom_point(data = MCMC_sum2, aes(log(avgbrood), mean.est.pmult, color = group), alpha = 0.2, size = 2.5) +
    geom_point(data = dat2, aes(log(avgbrood), pmult, color = group)) +
    #stat_smooth(data = dat, aes(log(avgbrood), pmult), color = "red") +
    labs(x="Log litter/clutch size", y="Probability of multiple paternity") +
    lims(y = c(0,1)) +
    scale_x_continuous(breaks = c(seq(0,15, by = 2)) , limits = c(0,15)) +
    theme_bw(base_size = 16)
)
dev.off()


# ################ DIAGNOSTICS BELOW ###################
# current <- 1 #change for different species (61 total)
# nam.current <- names(jags.results[current])
# jags.post <- jags.results[[nam.current]]
# gelman.diag(jags.post, multivariate = F) # values below 1.1 should be okay
# gelman.plot(jags.post, multivariate = F)
# n.eff <- effectiveSize(jags.post)
# 
# # visualize trace and posterior plots
# par(mar=c(2,2,1,1))
# plot(jags.post)
# 
# ##### make inference #####
# post.summ(jags.post, "q") # probability of success for each sire for a given data-gen run
# post.summ(jags.post, "nsires")
# post.summ(jags.post, "nmates")
# post.summ(jags.post, "avg.R") # avg number of sires given the est prob of success for each sire, using avg litter size
# post.summ(jags.post, "avg.M")
# post.summ(jags.post, "est.pmult.kq")

# #### Compare the est.sires (calculated with est. q) to true.sires
# sires.results <- sires.results[order(sires.results$species),]
# head(sires.results)
# for(i in unique(sires.results$true.sires)){
#   print(apply(sires.results[sires.results$true.sires==i, 2:5], 2, mean))
# }


################ CALCULATE RESIDUALS ###################
# MCMC_resids <- MCMC_summary
# loess_pmult <- loess(mean.est.pmult ~ avgbrood, data = MCMC_resids)
# MCMC_resids$pred_pmult <- predict(loess_pmult, newdata = MCMC_resids$avgbrood)
# # lower limit
# loess_pmult_lcl <- loess(est.pmult2.5 ~ avgbrood, data = MCMC_resids)
# MCMC_resids$pred_pmult_lcl <- predict(loess_pmult_lcl, newdata = MCMC_resids$avgbrood)
# # upper limit
# loess_pmult_ucl <- loess(est.pmult97.5 ~ avgbrood, data = MCMC_resids)
# MCMC_resids$pred_pmult_ucl <- predict(loess_pmult_ucl, newdata = MCMC_resids$avgbrood)
# # residual calculations
# MCMC_resids$resid_pmult <- MCMC_resids$pred_pmult - MCMC_resids$pmult
# MCMC_resids$resid_pmult_lcl <- MCMC_resids$pred_pmult_lcl - MCMC_resids$pmult
# MCMC_resids$resid_pmult_ucl <- MCMC_resids$pred_pmult_ucl - MCMC_resids$pmult
# 
# # save dataframe with all residual calculations
# save(MCMC_resids, file = "MCMC_resids.rda")
# write.csv(MCMC_resids, file = "MCMC_resids.csv")

## log residuals
MCMC_resids <- MCMC_summary
loess_pmult <- loess(mean.est.pmult ~ log(avgbrood), data = MCMC_resids)
MCMC_resids$pred_pmult <- predict(loess_pmult, newdata = log(MCMC_resids$avgbrood))
# lower limit
loess_pmult_lcl <- loess(est.pmult2.5 ~ log(avgbrood), data = MCMC_resids)
MCMC_resids$pred_pmult_lcl <- predict(loess_pmult_lcl, newdata = log(MCMC_resids$avgbrood))
# upper limit
loess_pmult_ucl <- loess(est.pmult97.5 ~ log(avgbrood), data = MCMC_resids)
MCMC_resids$pred_pmult_ucl <- predict(loess_pmult_ucl, newdata = log(MCMC_resids$avgbrood))
# residual calculations
MCMC_resids$resid_pmult <- MCMC_resids$pred_pmult - MCMC_resids$pmult
MCMC_resids$resid_pmult_lcl <- MCMC_resids$pred_pmult_lcl - MCMC_resids$pmult
MCMC_resids$resid_pmult_ucl <- MCMC_resids$pred_pmult_ucl - MCMC_resids$pmult

MCMC_resids$group <- c(rep("mammals", 49), rep("fish", 28), rep("herps", 16), rep("inverts", 30))
dat2 <- dat
dat2$group <- c(rep("mammals", 49), rep("fish", 28), rep("herps", 18), rep("inverts", 30))
dat25 <- dat2[!dat2$avgbrood>25000,]


save(MCMC_resids, file = "/Users/Hannah/Dropbox/Paternity2/BayesianMCMC/all_species_tog/MCMC_resids.rda")


## subset data to only include species with avgbrood <= 25k
MCMC_resids_25k <- MCMC_resids[!MCMC_resids$avgbrood>25000,]

## original version (not grouped)
pdf(file = paste0("p_CL.pdf"), width = 11, height = 7.5)
print(
  ggplot() +
    stat_smooth(data = MCMC_resids_25k, aes(log(avgbrood), pred_pmult), color = "red", se = FALSE, method = "loess") +
    geom_point(data = MCMC_resids_25k, aes(log(avgbrood), mean.est.pmult), color = "red", alpha = 0.5, size = 2.5) +
    stat_smooth(data = MCMC_resids_25k, aes(log(avgbrood), pred_pmult_ucl), color = "blue", se = FALSE, size = 0.5, method = "loess") +
    stat_smooth(data = MCMC_resids_25k, aes(log(avgbrood), pred_pmult_lcl), color = "blue", se = FALSE, size = 0.5, method = "loess") +
    stat_smooth(data = MCMC_resids_25k, aes(log(avgbrood), pred_pmult), color = "green3", se = FALSE, method = "lm") +
    geom_point(data = dat25, aes(log(avgbrood), pmult), alpha = 0.5, size = 2.5) +
    stat_smooth(data = dat25, aes(log(avgbrood), pmult), se = FALSE, color = "black", method = "lm") +
    labs(x="Log litter/clutch size", y="Probability of multiple paternity") +
    lims(y = c(0,1)) +
    scale_x_continuous(breaks = c(seq(0,11, by = 2)) , limits = c(0,11)) +
    theme_bw(base_size = 16)
)
dev.off()


# plot pred_pmult with lower and upper limits
mypal_jco <- pal_jco()(10)[c(8,9,7,1)]
library(RColorBrewer)
mypal <- brewer.pal(n=8, "Dark2")[c(1,2,4,3)]
library(wesanderson)
library(ggsci)
gg2 <- ggplot() +
  geom_point(data = MCMC_resids_25k, aes(log(avgbrood), mean.est.pmult, fill = group, shape = "Estimated"), alpha = 0.5, size = 2.5) +
  geom_point(data = dat25, aes(log(avgbrood), pmult, fill = group, shape = "Observed"), alpha = 0.7, size = 2.5) +
  # stat_smooth(data = dat, aes(log(avgbrood), pmult), se = FALSE, color = "black", method = "lm") +
  # stat_smooth(data = MCMC_resids_25k, aes(log(avgbrood), pred_pmult), color = "green3", se = FALSE, method = "lm") +
  stat_smooth(data = MCMC_resids_25k, aes(log(avgbrood), pred_pmult), color = "red", se = FALSE, method = "loess") +
  stat_smooth(data = MCMC_resids_25k, aes(log(avgbrood), pred_pmult_ucl), color = "blue", se = FALSE, size = 0.5, method = "loess") +
  stat_smooth(data = MCMC_resids_25k, aes(log(avgbrood), pred_pmult_lcl), color = "blue", se = FALSE, size = 0.5, method = "loess") +
  labs(x="Log litter/clutch size", y="Probability of multiple paternity") +
  lims(y = c(0,1)) +
  # scale_color_jco(name = "Taxon group") + 
  # scale_color_manual(values = wes_palette("Darjeeling1")[c(5,2,3,1)]) + 
  # scale_color_brewer(palette = "Dark2", type = "qual") +
  scale_fill_manual(name = "Taxa", values = mypal,
                    labels = c("Fish", "Reptiles & amphibians","Invertebrates", "Mammals")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  scale_shape_manual(name = "Type", labels = c("Estimated", "Observed"), 
                     values = c(21,23)) +
  scale_x_continuous(breaks = c(seq(0,11, by = 2)) , limits = c(0,11)) +
  theme_bw(base_size = 16)
pdf(file = paste0("p_CL_group.pdf"), width = 10, height = 6)
print(gg2)
dev.off()

ggsave("p_CL_group.eps", plot = gg2, device = cairo_ps, 
       width = 10, height = 6, dpi = 320)




############### GRAPHICAL ABSTRACT ###################
ga1 <- ggplot() +
  # geom_point(data = MCMC_resids_25k, aes(log(avgbrood), mean.est.pmult, fill = group, shape = "Estimated"), alpha = 0.5, size = 2.5) +
  geom_point(data = dat25, aes(log(avgbrood), pmult, fill = group), alpha = 0.7, size = 3.5, shape = 23) +
  # stat_smooth(data = dat, aes(log(avgbrood), pmult), se = FALSE, color = "black", method = "lm") +
  # stat_smooth(data = MCMC_resids_25k, aes(log(avgbrood), pred_pmult), color = "green3", se = FALSE, method = "lm") +
  geom_ribbon(data = MCMC_resids_25k, aes(x = log(avgbrood), ymin = pred_pmult_lcl, ymax = pred_pmult_ucl), fill = "grey60", alpha = 0.6) +
  stat_smooth(data = MCMC_resids_25k, aes(log(avgbrood), pred_pmult), color = "black", se = FALSE, method = "loess") +
  # stat_smooth(data = MCMC_resids_25k, aes(log(avgbrood), pred_pmult_ucl), color = "blue", se = FALSE, size = 0.5, method = "loess") +
  # stat_smooth(data = MCMC_resids_25k, aes(log(avgbrood), pred_pmult_lcl), color = "blue", se = FALSE, size = 0.5, method = "loess") +
  labs(x="Log litter/clutch size", y="Probability of multiple paternity") +
  lims(y = c(0,1)) +
  # scale_color_jco(name = "Taxon group") + 
  # scale_color_manual(values = wes_palette("Darjeeling1")[c(5,2,3,1)]) + 
  # scale_color_brewer(palette = "Dark2", type = "qual") +
  scale_fill_manual(name = "Taxa:", values = mypal, 
                    labels = c("Fish", "Reptiles & amphibians","Invertebrates", "Mammals")) +
  guides(fill = guide_legend(override.aes = list(shape = 23))) +
  # scale_shape_manual(name = "Type", labels = c("Estimated", "Observed"), 
  #                    values = c(21,23)) +
  scale_x_continuous(breaks = c(seq(0,11, by = 2)) , limits = c(0,11)) +
  theme_bw(base_size = 16) + 
  theme(legend.position="bottom", legend.box = "vertical", 
        legend.background = element_rect(size=0.5, linetype="solid", colour = "grey50"))
ggsave("GA2.eps", plot = ga1, device = cairo_ps, 
       width = 8, height = 7, dpi = 600)



############### GRAPHICAL ABSTRACT - OPTION 2 ############### 
ga2 <- ggplot(data = MCMC_resids_25k) +
  geom_segment(aes(x = log(avgbrood), xend = log(avgbrood), y = pred_pmult, yend = pmult, color = group),
               lty = "solid", size = 3, alpha = 0.2) +
  # geom_point(data = MCMC_resids_25k, aes(log(avgbrood), mean.est.pmult, fill = group, shape = "Estimated"), alpha = 0.5, size = 2.5) +
  geom_point(aes(log(avgbrood), pmult, fill = group), alpha = 0.7, size = 3.5, shape = 23) +
  geom_ribbon(aes(x = log(avgbrood), ymin = pred_pmult_lcl, ymax = pred_pmult_ucl), fill = "grey60", alpha = 0.6) +
  stat_smooth(aes(log(avgbrood), pred_pmult), color = "black", se = FALSE, method = "loess") +
  # stat_smooth(data = MCMC_resids_25k, aes(log(avgbrood), pred_pmult_ucl), color = "blue", se = FALSE, size = 0.5, method = "loess") +
  # stat_smooth(data = MCMC_resids_25k, aes(log(avgbrood), pred_pmult_lcl), color = "blue", se = FALSE, size = 0.5, method = "loess") +
  labs(x="Log litter/clutch size", y="Probability of multiple paternity") +
  lims(y = c(0,1)) +
  # scale_color_jco(name = "Taxon group") + 
  # scale_color_manual(values = wes_palette("Darjeeling1")[c(5,2,3,1)]) + 
  # scale_color_brewer(palette = "Dark2", type = "qual") +
  scale_fill_manual(name = "Taxa:", values = mypal, 
                    labels = c("Fish", "Reptiles & amphibians","Invertebrates", "Mammals")) +
  guides(fill = guide_legend(override.aes = list(shape = 23))) +
  scale_color_manual(name = "Taxa:", values = mypal, 
                    labels = c("Fish", "Reptiles & amphibians","Invertebrates", "Mammals")) +
  # scale_shape_manual(name = "Type", labels = c("Estimated", "Observed"), 
  #                    values = c(21,23)) +
  scale_x_continuous(breaks = c(seq(0,11, by = 2)) , limits = c(0,11)) +
  theme_bw(base_size = 16) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position="bottom", legend.box = "vertical", 
        legend.background = element_rect(size=0.5, linetype="solid", colour = "grey50"))


ggsave("GA2b.eps", plot = ga2, device = cairo_ps, 
       width = 8, height = 7, dpi = 600)



################ CALCULATE EFFECT SIZES ###################
# mamm <- read.csv("/Users/Hannah/Dropbox/Paternity2/BayesianMCMC/mammals_Avise/single_run/paternity_mammals_Avise.csv")
# fish <- read.csv("/Users/Hannah/Dropbox/Paternity2/BayesianMCMC/fish/single_run/paternity_fish.csv")
# herp <- read.csv("/Users/Hannah/Dropbox/Paternity2/BayesianMCMC/herps/single_run/paternity_herps.csv")
# invert <- read.csv("/Users/Hannah/Dropbox/Paternity2/BayesianMCMC/inverts/single_run/paternity_inverts.csv")
# dat <- rbind(mamm, fish, herp, invert)
# dat2 <- dat
# dat2$group <- c(rep("mammals", 49), rep("fish", 28), rep("herps", 18), rep("inverts", 30))
# dat25 <- dat2[!dat2$avgbrood>25000,]

dat_1 <- dat25[!is.na(dat25$avgsire),]

# load("/Users/Hannah/Dropbox/Paternity2/BayesianMCMC/all_species_tog/MCMC_resids.rda")
# MCMC_resids_25k <- MCMC_resids[!MCMC_resids$avgbrood>25000,]

# yi = Bayes pmult - true pmult (i.e. MCMC_total_resids$resid_pmult)
# vi = pmult*(1-pmult)/nbrood

dat_singlerun <- merge(dat_1[,c(1:3)], MCMC_resids_25k, by = c("species", "avgbrood"))
# order by brood size
dat_singlerun <- dat_singlerun[order(dat_singlerun$avgbrood, dat_singlerun$pmult, dat_singlerun$nbrood),]

### variance p(1-p)/n for fixed p assumption
dat_singlerun$vi <- dat_singlerun$pmult*(1-dat_singlerun$pmult)/dat_singlerun$nbrood

library(metafor)
paternity.meta.bayes <- rma(resid_pmult, vi, data = dat_singlerun)
paternity.meta.bayes
pdf(file = paste0("forest_bayes_singlerun-w-labels.pdf"), width = 12, height = 25)
print(
  forest(paternity.meta.bayes, slab = dat_singlerun$species, 
         order = order(dat_singlerun$avgbrood), cex = 1,
         xlim = c(-2.5,2),
         xlab = expression(paste(p[B], " - p"))),
  quote = FALSE,
  text(-2.5, 126, "Species",  pos=4),
  text(2, 126, expression(paste(p[B], " - p [95% CI]")), pos=2)
)
dev.off()


# k as moderator 
paternity.meta.k <- rma(resid_pmult, vi, mods = ~ avgbrood, data = dat_singlerun)
paternity.meta.k
pdf(file = paste0("forest_bayes_singlerun_k-mod.pdf"), width = 15, height = 25)
print(
  forest(paternity.meta.k, slab = dat_singlerun$species, 
         order = order(dat_singlerun$avgbrood), cex = 1)
)
dev.off()

# k as moderator (data without pmult = 0 or 1)
paternity.meta.k <- rma(resid_pmult, vi, mods = ~ avgbrood, data = dat_singlerun[!dat_singlerun$pmult==0 & !dat_singlerun$pmult==1,])
paternity.meta.k
forest(paternity.meta.k, slab = dat_singlerun[!dat_singlerun$pmult==0 & !dat_singlerun$pmult==1,]$species, 
       order = order(dat_singlerun[!dat_singlerun$pmult==0 & !dat_singlerun$pmult==1,]$avgbrood), cex = 1)



#### Forest plot with Hedge's g by animal group
datg <- within(dat_singlerun,
               {m1i <- mean.est.pmult
               m2i <- pmult
               sd1i <- sqrt(pmult*(1-pmult))
               n1i <- nbrood})

## fit overall random-effects model
paternity.meta.all <- rma(measure="SMD", 
                          m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd1i, n2i=n1i, 
                          data = datg[!datg$pmult==0 & !datg$pmult==1,])
paternity.meta.all


### fit meta-regression model to test for subgroup differences
paternity.meta.g <- rma(measure="SMD", 
                        m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd1i, n2i=n1i, 
                        mods = ~ factor(group), 
                        data = datg[!datg$pmult==0 & !datg$pmult==1,])
paternity.meta.g


### fit random-effects model in the four subgroups
res.m <- rma(measure="SMD", 
             m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd1i, n2i=n1i, 
             subset = (group=="mammals"), 
             data = datg[!datg$pmult==0 & !datg$pmult==1,])
res.h <- rma(measure="SMD", 
             m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd1i, n2i=n1i, 
             subset = (group=="herps"), 
             data = datg[!datg$pmult==0 & !datg$pmult==1,])
res.f <- rma(measure="SMD", 
             m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd1i, n2i=n1i, 
             subset = (group=="fish"), 
             data = datg[!datg$pmult==0 & !datg$pmult==1,])
res.i <- rma(measure="SMD", 
             m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd1i, n2i=n1i, 
             subset = (group=="inverts"), 
             data = datg[!datg$pmult==0 & !datg$pmult==1,])


## source custom forest.rma2(...) function for bold headers
source("forest.rma2.R")

cex1 <- 0.55

pdf(file = paste0("forest_all_singlerun.pdf"), width = 8, height = 15)
print(
  
  ## set up main forest plot structure
  forest.rma2(paternity.meta.all, xlim=c(-10, 11), cex=cex1,
              slab = datg[!datg$pmult==0 & !datg$pmult==1,]$species, 
              order = order(factor(datg[!datg$pmult==0 & !datg$pmult==1,]$group, 
                                   levels = c("inverts","fish","herps","mammals")), 
                            rev(datg[!datg$pmult==0 & !datg$pmult==1,]$avgbrood)), 
              rows=c(3:20,25:44,49:61,66:111), ylim=c(-1, 115),
              # xlab = expression(paste(p[B], " - p")), mlab = "",
              # header = c("Species", expression(paste(p[B] - p, "  [95% CI]"))),
              xlab = "Hedge's g  [95% CI]", mlab = "",
              header = c("Species", "Hedge's g  [95% CI]")),
  
  ### add text with Q-value, dfs, p-value, and I^2 statistic
  text(-10, -1, pos=4, cex=cex1, bquote(paste("RE Model for All Species (Q = ",
                                              .(formatC(paternity.meta.all$QE, digits=2, format="f")), ", df = ", .(paternity.meta.all$k - paternity.meta.all$p),
                                              ", p = ", .(formatC(paternity.meta.all$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                              .(formatC(paternity.meta.all$I2, digits=1, format="f")), "%)"))),
  
  ### add text for the test of subgroup differences
  text(-10, -2.2, pos=4, cex=cex1, bquote(paste("Test for Subgroup Differences: ",
                                                Q[M], " = ", .(formatC(paternity.meta.g$QM, digits=2, format="f")), ", df = ", .(paternity.meta.g$p - 1),
                                                ", p = ", .(formatC(paternity.meta.g$QMp, digits=2, format="f"))))),
  
  ### set font expansion factor (as in forest() above) and use bold italic
  ### font and save original settings in object 'op'
  op <- par(cex=cex1, font=4),
  
  ### add text for the subgroups
  text(-10, c(112,62,45,21), pos=4, c("Mammals", "Herps", "Fish", "Invertebrates")),
  
  ### switch to bold font
  # par(font=2),
  
  ## add headings text manually for expressions
  # text(c(-12,10), 116, pos=4, c("Species", expression(paste(p[B] - p, "  [95% CI]"))),
  
  ### set par back to the original settings
  par(op),
  
  ### add summary polygons for the four subgroups
  addpoly(res.m, row=64.5, cex=cex1, mlab=""),
  addpoly(res.h, row=47.5, cex=cex1, mlab=""),
  addpoly(res.f, row=23.5, cex=cex1, mlab=""),
  addpoly(res.i, row= 1.5, cex=cex1, mlab=""),
  
  ### add text with Q-value, dfs, p-value, and I^2 statistic for subgroups
  text(-10, 64.5, pos=4, cex=cex1, bquote(paste("RE Model for Mammals (Q = ",
                                                .(formatC(res.m$QE, digits=2, format="f")), ", df = ", .(res.m$k - res.m$p),
                                                ", p = ", .(formatC(res.m$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                                .(formatC(res.m$I2, digits=1, format="f")), "%)"))),
  text(-10, 47.5, pos=4, cex=cex1, bquote(paste("RE Model for Herps (Q = ",
                                                .(formatC(res.h$QE, digits=2, format="f")), ", df = ", .(res.h$k - res.h$p),
                                                ", p = ", .(formatC(res.h$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                                .(formatC(res.h$I2, digits=1, format="f")), "%)"))),
  text(-10, 23.5, pos=4, cex=cex1, bquote(paste("RE Model for Fish (Q = ",
                                                .(formatC(res.f$QE, digits=2, format="f")), ", df = ", .(res.f$k - res.f$p),
                                                ", p = ", .(formatC(res.f$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                                .(formatC(res.f$I2, digits=1, format="f")), "%)"))),
  text(-10, 1.5, pos=4, cex=cex1, font = 2, bquote(paste("RE Model for Inverts (Q = ",
                                               .(formatC(res.i$QE, digits=2, format="f")), ", df = ", .(res.i$k - res.i$p),
                                               ", p = ", .(formatC(res.i$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                               .(formatC(res.i$I2, digits=1, format="f")), "%)"))),
  # par(op),
  quote = FALSE, print.gap = NULL, right = FALSE, max = NULL, width = NULL, useSource = TRUE
)
dev.off()




################ E(M) AND E(S) VS K PLOT ###################
# ## Plot estimated M from data using Bayesian p overlaid on observed avg M vs obs. littersize
# loess_R <- loess(mean.avgR ~ avgbrood, data = MCMC_resids_25k)
# MCMC_resids_25k$pred_R <- predict(loess_R, newdata = MCMC_resids_25k$avgbrood)
# # lower limit
# loess_R_lcl <- loess(avgR2.5 ~ avgbrood, data = MCMC_resids_25k)
# MCMC_resids_25k$pred_R_lcl <- predict(loess_R_lcl, newdata = MCMC_resids_25k$avgbrood)
# # upper limit
# loess_R_ucl <- loess(avgR97.5 ~ avgbrood, data = MCMC_resids_25k)
# MCMC_resids_25k$pred_R_ucl <- predict(loess_R_ucl, newdata = MCMC_resids_25k$avgbrood)

## log residuals
loess_R <- loess(mean.avgR ~ log(avgbrood), data = MCMC_resids_25k)
MCMC_resids_25k$pred_R <- predict(loess_R, newdata = log(MCMC_resids_25k$avgbrood))
# lower limit
loess_R_lcl <- loess(avgR2.5 ~ log(avgbrood), data = MCMC_resids_25k)
MCMC_resids_25k$pred_R_lcl <- predict(loess_R_lcl, newdata = log(MCMC_resids_25k$avgbrood))
# upper limit
loess_R_ucl <- loess(avgR97.5 ~ log(avgbrood), data = MCMC_resids_25k)
MCMC_resids_25k$pred_R_ucl <- predict(loess_R_ucl, newdata = log(MCMC_resids_25k$avgbrood))

# MCMC_resids_25k$group <- c(rep("mammals", 49), rep("fish", 28), rep("herps", 16), rep("inverts", 30))
# dat2 <- dat
# dat2$group <- c(rep("mammals", 49), rep("fish", 28), rep("herps", 18), rep("inverts", 30))


pdf(file = paste0("R_CL.pdf"), width = 10, height = 6)
print(
  ggplot() +
    geom_point(data = MCMC_resids_25k, aes(log(avgbrood), mean.avgR, fill = group, shape = "Estimated"), alpha = 0.5, size = 2.5) +
    geom_point(data = dat25, aes(log(avgbrood), avgsire, fill = group, shape = "Observed"), alpha = 0.5, size = 2.5) +
    stat_smooth(data = MCMC_resids_25k, aes(log(avgbrood), pred_R), color = "red", se = FALSE, method = "loess") +
    stat_smooth(data = MCMC_resids_25k, aes(log(avgbrood), pred_R_ucl), color = "blue", se = FALSE, size = 0.5, method = "loess") +
    stat_smooth(data = MCMC_resids_25k, aes(log(avgbrood), pred_R_lcl), color = "blue", se = FALSE, size = 0.5, method = "loess") +
    # stat_smooth(data = dat25, aes(log(avgbrood), avgsire), color = "black", method = "lm", se = FALSE) +
    scale_fill_manual(name = "Taxa", values = mypal,
                      labels = c("Fish", "Reptiles & amphibians","Invertebrates", "Mammals")) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    scale_shape_manual(name = "Type", labels = c("Estimated", "Observed"), 
                       values = c(21,23)) +
    labs(x="Log litter/clutch size", y="Number of sires") +
    scale_x_continuous(breaks = c(seq(0,11, by = 2)), limits = c(0,11)) +
    scale_y_continuous(breaks = c(seq(0,20, by = 5)), limits = c(0,20)) +
    theme_bw(base_size = 16)
)
dev.off()


## Replicate Ash's plot of #mates/sires vs brood size with Bayesian estimates
dat_1 <- dat25[!is.na(dat25$avgsire),]
# Combinatorics avgmates calculated using pred_pmult and true k
dat_1$avgmate <- exp(log(1-dat_1$pmult)/(1-dat_1$avgbrood))

# Used the Avise data avgbrood and avgsire to obtain avgR = k*q/(1-(1-q)^k) and avgM = r/q
# calculate predicted M
loess_M <- loess(mean.avgM ~ log(avgbrood), data = MCMC_resids_25k)
MCMC_resids_25k$pred_M <- predict(loess_M, newdata = log(MCMC_resids_25k$avgbrood))
# plot
pdf(file = paste0("Avise_MS-k.pdf"), width = 12, height = 7.5)
print(
  ggplot() +
    # avgR in red is calculated from mean of ZTB dist, since S ~ ZTB
    stat_smooth(data = MCMC_resids_25k, aes(log(avgbrood), log(pred_R), linetype = "Estimated # sires using\ndata K and Bayesian q"), color = "red", method = "loess", se = TRUE) +
    # points from which the above line was generated:
    #geom_point(data = MCMC_resids_25k, aes(avgbrood, mean.avgR), color = "red", alpha = 0.2) +
    # True number of average sires from the Avise data
    geom_point(data = dat_1, aes(log(avgbrood), log(avgsire), color = "Observed number of sires"), alpha = 0.5, size = 2.5) +
    # Avise used avgsires as estimated mates (i.e. he does not distinguish b/w mates and sires?)
    stat_smooth(data = dat_1, aes(log(avgbrood), log(avgsire), linetype = "Observed number of sires"), color = "black", method = "lm", se = FALSE) +
    # avgM in blue is calculated from mean of NB dist knowing S and q, since M ~ NB
    stat_smooth(data = MCMC_resids_25k, aes(log(avgbrood), log(pred_M), linetype = "Estimated # mates using\ndata S and Bayesian q"), color = "blue", method = "lm", se = TRUE) +
    # points from which the above line was generated:
    #geom_point(data = MCMC_resids_25k, aes(avgbrood, mean.avgM), color = "blue", alpha = 0.2) +
    # could instead plot avgmates calculated using true pmult and true k from combinatorics formula?
    stat_smooth(data = dat_1[!is.infinite(dat_1$avgmate),], aes(log(avgbrood), log(avgmate), linetype = "Estimated # mates using\ncombinatorial formula"), 
                color = "green3", method = "lm", se = TRUE) +
    geom_point(data = dat_1[!is.infinite(dat_1$avgmate),], aes(log(avgbrood), log(avgmate), color = "Estimated # mates using\ncombinatorial formula"), alpha = 0.5, size = 2.5) +
    labs(x="Log litter/clutch size", y="Log number of mates/sires") +
    scale_x_continuous(breaks = c(seq(0,11, by = 1)) , limits = c(0,11)) +
    scale_y_continuous(breaks = c(seq(0,16, by = 2)) , limits = c(-0.5,16)) +
    scale_colour_manual(name = "", values = c("Observed number of sires" = "black",
                                              "Estimated # mates using\ncombinatorial formula" = "green3")) +
    scale_linetype_manual(values = c("Estimated # mates using\ndata S and Bayesian q" = 1, # blue line
                                     "Estimated # mates using\ncombinatorial formula" = 1, # dark green line
                                     "Estimated # sires using\ndata K and Bayesian q" = 1, # red line
                                     "Observed number of sires" = 1), # black line
                          name = "",
                          guide = guide_legend(override.aes = list(color = c("green3", "blue", "red", "black"), 
                                                                   fill = c(NA,NA,NA,NA)))) +
    theme_bw(base_size = 16) + theme(legend.key.height=unit(15, "mm"))
)
dev.off()


# Used the generated data to create nsires and calculate estM = nsires/q
pdf(file = paste0("Bayes_gen_MS-k.pdf"), width = 12, height = 7.5)
print(
  ggplot() +
    # NB estimated mates using different values of q
    stat_function(fun = function(x) x/(1-0.2^x), aes(x = x, linetype = "Expected # mates from NB distribution"), data = data.frame(x = c(2,10)), color = "magenta", size = 1.2) +
    stat_function(fun = function(x) x/(1-0.3^x), aes(x = x, linetype = "Expected # mates from NB distribution"), data = data.frame(x = c(2,10)), color = "magenta", size = 1.2) +
    stat_function(fun = function(x) x/(1-0.4^x), aes(x = x, linetype = "Expected # mates from NB distribution"), data = data.frame(x = c(2,10)), color = "magenta", size = 1.2) +
    stat_function(fun = function(x) x/(1-0.5^x), aes(x = x, linetype = "Expected # mates from NB distribution"), data = data.frame(x = c(2,10)), color = "magenta", size = 1.2) +
    # # nsires in red is generated from Poission dist with lambda = avgsire in Avise data
    # stat_smooth(data = MCMC_resids_25k, aes(avgbrood, mean.nsires, linetype = "Generated # sires using Poisson(avgsire)"), color = "red", method = "lm", se = TRUE) +
    # # points from which the above line was generated:
    # geom_point(data = MCMC_resids_25k, aes(avgbrood, mean.nsires, color = "Generated # sires using Poisson(avgsire)"), alpha = 0.5) +
    # True number of average sires from the Avise data
    geom_point(data = dat_1, aes(log(avgbrood), log(avgsire), color = "Observed number of sires"), alpha = 0.5, size = 2.5) +
    # Avise used avgsires as estimated mates (i.e. he does not distinguish b/w mates and sires?)
    stat_smooth(data = dat_1, aes(log(avgbrood), log(avgsire), linetype = "Observed number of sires"), color = "black", method = "lm", se = FALSE) +
    # estM in blue is calculated from mean of NB dist knowing generated S and q, since M ~ NB
    stat_smooth(data = MCMC_resids_25k, aes(log(avgbrood), log(mean.estM), linetype = "Bayesian MCMC estimated # mates"), color = "blue", method = "lm", se = FALSE) +
    # points from which the above line was generated:
    geom_point(data = MCMC_resids_25k, aes(log(avgbrood), log(mean.estM), color = "Bayesian MCMC estimated # mates"), alpha = 0.5, size = 2.5) +
    # could instead plot avgmates calculated using true pmult and true k from combinatorics formula?
    stat_smooth(data = dat_1[!is.infinite(dat_1$avgmate),], aes(log(avgbrood), log(avgmate), linetype = "Estimated # mates using\ncombinatorial formula"), 
                color = "green3", method = "lm", se = FALSE) +
    geom_point(data = dat_1[!is.infinite(dat_1$avgmate),], aes(log(avgbrood), log(avgmate), color = "Estimated # mates using\ncombinatorial formula"), alpha = 0.5, size = 2.5) +
    labs(x="Log litter/clutch size", y="Log number of mates/sires") +
    scale_x_continuous(breaks = c(seq(1,11, by = 1)) , limits = c(min(log(dat_1$avgbrood)), max(log(dat_1$avgbrood)))) +
    scale_y_continuous(breaks = c(seq(0,16, by = 2)) , limits = c(-0.5,16)) +
    scale_colour_manual(name = "", values = c("Bayesian MCMC estimated # mates" = "blue",
                                              "Estimated # mates using\ncombinatorial formula" = "green3",
                                              #"Generated # sires using Poisson(avgsire)" = "red",
                                              "Observed number of sires" = "black")) +
    scale_linetype_manual(values = c("Bayesian MCMC estimated # mates" = 1, # blue line,
                                     "Estimated # mates using\ncombinatorial formula" = 1, # green3 line
                                     "Expected # mates from NB distribution" = 1, # magenta line
                                     #"Generated # sires using Poisson(avgsire)" = 1, # red line
                                     "Observed number of sires" = 1), # black line
                          name = "",
                          guide = guide_legend(override.aes = list(color = c("blue","green3", 
                                                                             "magenta", 
                                                                             #"red", 
                                                                             "black"), 
                                                                   fill = c(NA,NA,
                                                                            NA,
                                                                            #NA,
                                                                            NA)))) +
    theme_bw(base_size = 16) + theme(legend.key.height=unit(10, "mm"))
)
dev.off()




##### Clustering - EXPERIMENTAL - NOT USED #####
MCMC_resids2 <- MCMC_resids
colnames(MCMC_resids2)
MCMC_resids2$group <- c(rep("mammals", 49), rep("fish", 28), rep("herps", 16), rep("inverts", 30))
MCMC_group1 <- MCMC_resids2[,c("avgbrood", "avgsire", "mean.q", "resid_pmult", "mean.estM")]
MCMC_group <- MCMC_group1[MCMC_group1$avgbrood<25000,]

# fit <- kmeans(MCMC_group, 5)
# # library(cluster)
# clusplot(MCMC_group, fit$cluster, color=TRUE, shade=TRUE, labels=2, lines=0)
# # library(fpc)
# plotcluster(MCMC_group, fit$cluster)


## Model-based clustering
library(mclust)
fit <- Mclust(MCMC_group)
plot(fit) # plot results 
summary(fit)


mod2 <- MclustDA(MCMC_resids2[,c("avgbrood", "avgsire", "mean.q", 
                                 "resid_pmult", "mean.estM", "group")], MCMC_resids2$group)
summary(mod2)
plot(mod2)



#### Compare correlations
library(cocor)

## Mammals
mamdat <- MCMC_resids2[MCMC_resids2$group=="mammals",]

## Avies's orginal result
cor.test(mamdat$avgbrood, mamdat$pmult)

## Testing Avise's correlation against our "null correlation"
cocor(~ pmult + avgbrood | mean.est.pmult + avgbrood, data = mamdat, test = "steiger1980",
      alternative="two.sided", alpha=0.05, conf.level=0.95, null.value=0)

cocor.dep.groups.overlap(r.jk = cor(mamdat$avgbrood, mamdat$pmult, method = "spearman"), 
                         r.jh = cor(mamdat$avgbrood, mamdat$mean.est.pmult, method = "spearman"), 
                         r.kh = cor(mamdat$pmult, mamdat$mean.est.pmult, method = "spearman"), 
                         n = nrow(mamdat), alternative = "two.sided",
                         test = "steiger1980", alpha = 0.05, conf.level = 0.95, null.value = 0,
                         data.name = NULL, var.labels = NULL, return.htest = FALSE)


## Herps
herpdat <- MCMC_resids2[MCMC_resids2$group=="herps",]

## Avies's orginal result
cor.test(herpdat$avgbrood, herpdat$pmult)

## Testing Avise's correlation against our "null correlation"
cocor(~ pmult + avgbrood | mean.est.pmult + avgbrood, data = herpdat, test = "steiger1980",
      alternative="two.sided", alpha=0.05, conf.level=0.95, null.value=0)

cocor.dep.groups.overlap(r.jk = cor(herpdat$avgbrood, herpdat$pmult, method = "spearman"), 
                         r.jh = cor(herpdat$avgbrood, herpdat$mean.est.pmult, method = "spearman"), 
                         r.kh = cor(herpdat$pmult, herpdat$mean.est.pmult, method = "spearman"), 
                         n = nrow(herpdat), alternative = "two.sided",
                         test = "steiger1980", alpha = 0.05, conf.level = 0.95, null.value = 0,
                         data.name = NULL, var.labels = NULL, return.htest = FALSE)


## Fish
fishdat <- MCMC_resids2[MCMC_resids2$group=="fish" & MCMC_resids2$avgbrood<25000,]

## Avies's orginal result
cor.test(fishdat$avgbrood, fishdat$pmult)

## Testing Avise's correlation against our "null correlation"
cocor(~ pmult + avgbrood | mean.est.pmult + avgbrood, data = fishdat, test = "steiger1980",
      alternative="two.sided", alpha=0.05, conf.level=0.95, null.value=0)

cocor.dep.groups.overlap(r.jk = cor(fishdat$avgbrood, fishdat$pmult, method = "spearman"), 
                         r.jh = cor(fishdat$avgbrood, fishdat$mean.est.pmult, method = "spearman"), 
                         r.kh = cor(fishdat$pmult, fishdat$mean.est.pmult, method = "spearman"), 
                         n = nrow(fishdat), alternative = "two.sided",
                         test = "steiger1980", alpha = 0.05, conf.level = 0.95, null.value = 0,
                         data.name = NULL, var.labels = NULL, return.htest = FALSE)


## Inverts
invdat <- MCMC_resids2[MCMC_resids2$group=="inverts" & MCMC_resids2$avgbrood<25000,]

## Avies's orginal result
cor.test(invdat$avgbrood, invdat$pmult)

## Testing Avise's correlation against our "null correlation"
cocor(~ pmult + avgbrood | mean.est.pmult + avgbrood, data = invdat, test = "steiger1980",
      alternative="two.sided", alpha=0.05, conf.level=0.95, null.value=0)

cocor.dep.groups.overlap(r.jk = cor(invdat$avgbrood, invdat$pmult, method = "spearman"), 
                         r.jh = cor(invdat$avgbrood, invdat$mean.est.pmult, method = "spearman"), 
                         r.kh = cor(invdat$pmult, invdat$mean.est.pmult, method = "spearman"), 
                         n = nrow(invdat), alternative = "two.sided",
                         test = "steiger1980", alpha = 0.05, conf.level = 0.95, null.value = 0,
                         data.name = NULL, var.labels = NULL, return.htest = FALSE)





#### Slope & correlation comparisons of Avise's results to our null ####
nulldat <- MCMC_resids2[, c("group", "avgbrood", "avgsire", "mean.est.pmult")]
colnames(nulldat)[4] <- "pmult"
nulldat$dummy <- "nullmod"
truedat <- MCMC_resids2[, c("group", "avgbrood", "avgsire", "pmult")]
truedat$dummy <- "truedat"

lmdat <- rbind(nulldat, truedat)

#### fit LM to null model and true data using dummy for each taxa group
## mammals
mfit <- lm(pmult ~ avgbrood*dummy, data = lmdat[lmdat$group=="mammals",])
summary(mfit)
anova(lm(pmult ~ avgbrood, data = lmdat[lmdat$group=="mammals",]), mfit)

## herps
hfit <- lm(pmult ~ avgbrood*dummy, data = lmdat[lmdat$group=="herps",])
summary(hfit)
anova(lm(pmult ~ avgbrood, data = lmdat[lmdat$group=="herps",]), hfit)

## fish
ffit <- lm(pmult ~ avgbrood*dummy, data = lmdat[lmdat$group=="fish",])
summary(ffit)
anova(lm(pmult ~ avgbrood, data = lmdat[lmdat$group=="fish",]), ffit)
## less than 25k offspring
ffit25k <- lm(pmult ~ avgbrood*dummy, data = lmdat[lmdat$group=="fish" & lmdat$avgbrood<25000,])
summary(ffit25k)
anova(lm(pmult ~ avgbrood, data = lmdat[lmdat$group=="fish" & lmdat$avgbrood<25000,]), ffit25k)

## inverts
ifit <- lm(pmult ~ avgbrood*dummy, data = lmdat[lmdat$group=="inverts",])
summary(ifit)
anova(lm(pmult ~ avgbrood, data = lmdat[lmdat$group=="inverts",]), ifit)
## less than 25k offspring
ifit25k <- lm(pmult ~ avgbrood*dummy, data = lmdat[lmdat$group=="inverts" & lmdat$avgbrood<25000,])
summary(ifit25k)
anova(lm(pmult ~ avgbrood, data = lmdat[lmdat$group=="inverts" & lmdat$avgbrood<25000,]), ifit25k)


library(ggplot2)
library(gridExtra)
## Mammals
pm <- ggplot(data = lmdat[lmdat$group=="mammals",], aes(x = avgbrood, y = pmult, group = dummy)) + 
  geom_point(aes(colour = dummy)) + 
  scale_color_manual(values = c("#00AFBB", "#FC4E07"), 
                     label = c("Null model", "True data"), name = "") +
  # xlim(0,400) + 
  geom_abline(intercept = mfit$coefficients[[1]], 
              slope = mfit$coefficients[[2]], color = "#00AFBB") + 
  geom_abline(intercept = (mfit$coefficients[[1]] + mfit$coefficients[[3]]), 
              slope = (mfit$coefficients[[2]] + mfit$coefficients[[4]]), 
              color = "#FC4E07") +
  theme_bw() + ggtitle(paste0("Mammals - Difference in slopes: p=", round(summary(mfit)$coefficients[4,4], 4)))


## Herps
ph <- ggplot(data = lmdat[lmdat$group=="herps",], aes(x = avgbrood, y = pmult, group = dummy)) + 
  geom_point(aes(colour = dummy)) + 
  scale_color_manual(values = c("#00AFBB", "#FC4E07"), 
                     label = c("Null model", "True data"), name = "") +
  # xlim(0,400) + 
  geom_abline(intercept = hfit$coefficients[[1]], 
              slope = hfit$coefficients[[2]], color = "#00AFBB") + 
  geom_abline(intercept = (hfit$coefficients[[1]] + hfit$coefficients[[3]]), 
              slope = (hfit$coefficients[[2]] + hfit$coefficients[[4]]), 
              color = "#FC4E07") +
  theme_bw() + ggtitle(paste0("Herps - Difference in slopes: p=", round(summary(hfit)$coefficients[4,4], 4)))


## Fish (all)
pf <- ggplot(data = lmdat[lmdat$group=="fish",], aes(x = avgbrood, y = pmult, group = dummy)) + 
  geom_point(aes(colour = dummy)) + 
  scale_color_manual(values = c("#00AFBB", "#FC4E07"), 
                     label = c("Null model", "True data"), name = "") +
  # xlim(0,400) + 
  geom_abline(intercept = ffit$coefficients[[1]], 
              slope = ffit$coefficients[[2]], color = "#00AFBB") + 
  geom_abline(intercept = (ffit$coefficients[[1]] + ffit$coefficients[[3]]), 
              slope = (ffit$coefficients[[2]] + ffit$coefficients[[4]]), 
              color = "#FC4E07") +
  theme_bw() + ggtitle(paste0("Fish (all) - Difference in slopes: p=", round(summary(ffit)$coefficients[4,4], 4)))

## Fish (less than 25k broodsize)
pf25 <- ggplot(data = lmdat[lmdat$group=="fish" & lmdat$avgbrood<25000,], 
               aes(x = avgbrood, y = pmult, group = dummy)) + 
  geom_point(aes(colour = dummy)) + 
  scale_color_manual(values = c("#00AFBB", "#FC4E07"), 
                     label = c("Null model", "True data"), name = "") +
  # xlim(0,400) + 
  geom_abline(intercept = ffit25k$coefficients[[1]], 
              slope = ffit25k$coefficients[[2]], color = "#00AFBB") + 
  geom_abline(intercept = (ffit25k$coefficients[[1]] + ffit25k$coefficients[[3]]), 
              slope = (ffit25k$coefficients[[2]] + ffit25k$coefficients[[4]]), 
              color = "#FC4E07") +
  theme_bw() + ggtitle(paste0("Fish (< 25k avgbrood) - Difference in slopes: p=", round(summary(ffit25k)$coefficients[4,4], 4)))


## Inverts
pi <- ggplot(data = lmdat[lmdat$group=="inverts",], aes(x = avgbrood, y = pmult, group = dummy)) + 
  geom_point(aes(colour = dummy)) + 
  scale_color_manual(values = c("#00AFBB", "#FC4E07"), 
                     label = c("Null model", "True data"), name = "") +
  # xlim(0,400) + 
  geom_abline(intercept = ifit$coefficients[[1]], 
              slope = ifit$coefficients[[2]], color = "#00AFBB") + 
  geom_abline(intercept = (ifit$coefficients[[1]] + ifit$coefficients[[3]]), 
              slope = (ifit$coefficients[[2]] + ifit$coefficients[[4]]), 
              color = "#FC4E07") +
  theme_bw() + ggtitle(paste0("Inverts - Difference in slopes: p=", round(summary(ifit)$coefficients[4,4], 4)))

## Inverts (less than 25k broodsize)
pi25 <- ggplot(data = lmdat[lmdat$group=="inverts" & lmdat$avgbrood<25000,], 
               aes(x = avgbrood, y = pmult, group = dummy)) + 
  geom_point(aes(colour = dummy)) + 
  scale_color_manual(values = c("#00AFBB", "#FC4E07"), 
                     label = c("Null model", "True data"), name = "") +
  # xlim(0,400) + 
  geom_abline(intercept = ifit25k$coefficients[[1]], 
              slope = ifit25k$coefficients[[2]], color = "#00AFBB") + 
  geom_abline(intercept = (ifit25k$coefficients[[1]] + ifit25k$coefficients[[3]]), 
              slope = (ifit25k$coefficients[[2]] + ifit25k$coefficients[[4]]), 
              color = "#FC4E07") +
  theme_bw() + ggtitle(paste0("Inverts (< 25k avgbrood) - Difference in slopes: p=", round(summary(ifit25k)$coefficients[4,4], 4)))




## All taxa together (less than 25k broodsize)

fit25k <- lm(pmult ~ avgbrood*dummy, data = lmdat[lmdat$avgbrood<25000,])
summary(fit25k)
anova(lm(pmult ~ avgbrood, data = lmdat[lmdat$avgbrood<25000,]), fit25k)

pall <- ggplot(data = lmdat[lmdat$avgbrood<25000,], 
               aes(x = avgbrood, y = pmult, group = dummy)) + 
  geom_point(aes(colour = dummy)) + 
  scale_color_manual(values = c("#00AFBB", "#FC4E07"), 
                     label = c("Null model", "True data"), name = "") +
  # xlim(0,400) + 
  geom_abline(intercept = fit25k$coefficients[[1]], 
              slope = fit25k$coefficients[[2]], color = "#00AFBB") + 
  geom_abline(intercept = (fit25k$coefficients[[1]] + fit25k$coefficients[[3]]), 
              slope = (fit25k$coefficients[[2]] + fit25k$coefficients[[4]]), 
              color = "#FC4E07") +
  theme_bw() + ggtitle(paste0("All taxa (< 25k avgbrood) - Difference in slopes: p=", round(summary(fit25k)$coefficients[4,4], 4)))


grid.arrange(pm, ph, pf25, pi25, pall, ncol = 2)



#### fit rank LM to null model and true data using dummy for each taxa group
library(Rfit)
## mammals
mfit <- rfit(pmult ~ avgbrood*dummy, data = lmdat[lmdat$group=="mammals",])
summary(mfit)

## herps
hfit <- rfit(pmult ~ avgbrood*dummy, data = lmdat[lmdat$group=="herps",])
summary(hfit)

## fish
ffit <- rfit(pmult ~ avgbrood*dummy, data = lmdat[lmdat$group=="fish",])
summary(ffit)
## less than 25k offspring
ffit25k <- rfit(pmult ~ avgbrood*dummy, data = lmdat[lmdat$group=="fish" & lmdat$avgbrood<25000,])
summary(ffit25k)

## inverts
ifit <- rfit(pmult ~ avgbrood*dummy, data = lmdat[lmdat$group=="inverts",])
summary(ifit)
## less than 25k offspring
ifit25k <- rfit(pmult ~ avgbrood*dummy, data = lmdat[lmdat$group=="inverts" & lmdat$avgbrood<25000,])
summary(ifit25k)


library(ggplot2)
## Mammals
pm <- ggplot(data = lmdat[lmdat$group=="mammals",], aes(x = avgbrood, y = pmult, group = dummy)) + 
  geom_point(aes(colour = dummy)) + 
  scale_color_manual(values = c("#00AFBB", "#FC4E07"), 
                     label = c("Null model", "True data"), name = "") +
  # xlim(0,400) + 
  geom_abline(intercept = mfit$coefficients[[1]], 
              slope = mfit$coefficients[[2]], color = "#00AFBB") + 
  geom_abline(intercept = (mfit$coefficients[[1]] + mfit$coefficients[[3]]), 
              slope = (mfit$coefficients[[2]] + mfit$coefficients[[4]]), 
              color = "#FC4E07") +
  theme_bw() + ggtitle(paste0("Mammals - Difference in slopes: p=", round(summary(mfit)$coefficients[4,4], 4)))


## Herps
ph <- ggplot(data = lmdat[lmdat$group=="herps",], aes(x = avgbrood, y = pmult, group = dummy)) + 
  geom_point(aes(colour = dummy)) + 
  scale_color_manual(values = c("#00AFBB", "#FC4E07"), 
                     label = c("Null model", "True data"), name = "") +
  # xlim(0,400) + 
  geom_abline(intercept = hfit$coefficients[[1]], 
              slope = hfit$coefficients[[2]], color = "#00AFBB") + 
  geom_abline(intercept = (hfit$coefficients[[1]] + hfit$coefficients[[3]]), 
              slope = (hfit$coefficients[[2]] + hfit$coefficients[[4]]), 
              color = "#FC4E07") +
  theme_bw() + ggtitle(paste0("Herps - Difference in slopes: p=", round(summary(hfit)$coefficients[4,4], 4)))


## Fish (all)
pf <- ggplot(data = lmdat[lmdat$group=="fish",], aes(x = avgbrood, y = pmult, group = dummy)) + 
  geom_point(aes(colour = dummy)) + 
  scale_color_manual(values = c("#00AFBB", "#FC4E07"), 
                     label = c("Null model", "True data"), name = "") +
  # xlim(0,400) + 
  geom_abline(intercept = ffit$coefficients[[1]], 
              slope = ffit$coefficients[[2]], color = "#00AFBB") + 
  geom_abline(intercept = (ffit$coefficients[[1]] + ffit$coefficients[[3]]), 
              slope = (ffit$coefficients[[2]] + ffit$coefficients[[4]]), 
              color = "#FC4E07") +
  theme_bw() + ggtitle(paste0("Fish (all) - Difference in slopes: p=", round(summary(ffit)$coefficients[4,4], 4)))

## Fish (less than 25k broodsize)
pf25 <- ggplot(data = lmdat[lmdat$group=="fish" & lmdat$avgbrood<25000,], 
               aes(x = avgbrood, y = pmult, group = dummy)) + 
  geom_point(aes(colour = dummy)) + 
  scale_color_manual(values = c("#00AFBB", "#FC4E07"), 
                     label = c("Null model", "True data"), name = "") +
  # xlim(0,400) + 
  geom_abline(intercept = ffit25k$coefficients[[1]], 
              slope = ffit25k$coefficients[[2]], color = "#00AFBB") + 
  geom_abline(intercept = (ffit25k$coefficients[[1]] + ffit25k$coefficients[[3]]), 
              slope = (ffit25k$coefficients[[2]] + ffit25k$coefficients[[4]]), 
              color = "#FC4E07") +
  theme_bw() + ggtitle(paste0("Fish (< 25k avgbrood) - Difference in slopes: p=", round(summary(ffit25k)$coefficients[4,4], 4)))


## Inverts
pi <- ggplot(data = lmdat[lmdat$group=="inverts",], aes(x = avgbrood, y = pmult, group = dummy)) + 
  geom_point(aes(colour = dummy)) + 
  scale_color_manual(values = c("#00AFBB", "#FC4E07"), 
                     label = c("Null model", "True data"), name = "") +
  # xlim(0,400) + 
  geom_abline(intercept = ifit$coefficients[[1]], 
              slope = ifit$coefficients[[2]], color = "#00AFBB") + 
  geom_abline(intercept = (ifit$coefficients[[1]] + ifit$coefficients[[3]]), 
              slope = (ifit$coefficients[[2]] + ifit$coefficients[[4]]), 
              color = "#FC4E07") +
  theme_bw() + ggtitle(paste0("Inverts - Difference in slopes: p=", round(summary(ifit)$coefficients[4,4], 4)))

## Inverts (less than 25k broodsize)
pi25 <- ggplot(data = lmdat[lmdat$group=="inverts" & lmdat$avgbrood<25000,], 
               aes(x = avgbrood, y = pmult, group = dummy)) + 
  geom_point(aes(colour = dummy)) + 
  scale_color_manual(values = c("#00AFBB", "#FC4E07"), 
                     label = c("Null model", "True data"), name = "") +
  # xlim(0,400) + 
  geom_abline(intercept = ifit25k$coefficients[[1]], 
              slope = ifit25k$coefficients[[2]], color = "#00AFBB") + 
  geom_abline(intercept = (ifit25k$coefficients[[1]] + ifit25k$coefficients[[3]]), 
              slope = (ifit25k$coefficients[[2]] + ifit25k$coefficients[[4]]), 
              color = "#FC4E07") +
  theme_bw() + ggtitle(paste0("Inverts (< 25k avgbrood) - Difference in slopes: p=", round(summary(ifit25k)$coefficients[4,4], 4)))




## All taxa together (less than 25k broodsize)
fit25k <- rfit(pmult ~ avgbrood*dummy, data = lmdat[lmdat$avgbrood<25000,])
summary(fit25k)

pall <- ggplot(data = lmdat[lmdat$avgbrood<25000,], 
               aes(x = avgbrood, y = pmult, group = dummy)) + 
  geom_point(aes(colour = dummy)) + 
  scale_color_manual(values = c("#00AFBB", "#FC4E07"), 
                     label = c("Null model", "True data"), name = "") +
  # xlim(0,400) + 
  geom_abline(intercept = fit25k$coefficients[[1]], 
              slope = fit25k$coefficients[[2]], color = "#00AFBB") + 
  geom_abline(intercept = (fit25k$coefficients[[1]] + fit25k$coefficients[[3]]), 
              slope = (fit25k$coefficients[[2]] + fit25k$coefficients[[4]]), 
              color = "#FC4E07") +
  theme_bw() + ggtitle(paste0("All taxa (< 25k avgbrood) - Difference in slopes: p=", round(summary(fit25k)$coefficients[4,4], 4)))



grid.arrange(pm, ph, pf25, pi25, pall, ncol = 2)


