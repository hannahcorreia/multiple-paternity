###################################################### 
# Written by: Hannah Correia, Auburn University, USA #
#             Asheber Abebe, Auburn University, USA  #
######################################################

# clear the workspace
rm(list = ls(all = T))

library(tidyverse)
library(rjags)
library(R2OpenBUGS)
library(coda)
library(extraDistr)
library(dplyr)
# library(ggplot2)

##### specify post.summ function
## Function written by Ben Staton, Columbia River Inter-Tribal Fish Commission, Portland, OR
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
mamm <- read.csv("paternity_mammals_Avise.csv")
fish <- read.csv("paternity_fish.csv")
herp <- read.csv("paternity_herps.csv")
invert <- read.csv("paternity_inverts.csv")
invert <- invert %>% rename(nbrood = nbroods)
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
    est.pmult.kq[j] <- 1 - avgbrood[j] * q * (1 - q)^(avgbrood[j] - 1) / (1 - (1 - q)^avgbrood[j])
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



################ CALCULATE RESIDUALS ###################

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


save(MCMC_resids, file = "MCMC_resids.rda")


## subset data to only include species with avgbrood <= 25k
MCMC_resids_25k <- MCMC_resids[!MCMC_resids$avgbrood>25000,]


# plot pred_pmult with lower and upper limits (grouped by taxa)
library(RColorBrewer)
mypal <- brewer.pal(n=8, "Dark2")[c(1,2,4,3)]

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



################ CALCULATE EFFECT SIZES ###################
dat_1 <- dat25[!is.na(dat25$avgsire),]

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



#### Forest plot with Hedges' g by animal group
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
              xlab = "Hedges' g  [95% CI]", mlab = "",
              header = c("Species", "Hedges' g  [95% CI]")),
  
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


