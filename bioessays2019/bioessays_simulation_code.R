###################################################### 
# Written by: Hannah Correia, Auburn University, USA #
#             Asheber Abebe, Auburn University, USA  #
######################################################

library(rjags)
library(R2OpenBUGS)
library(coda)
library(extraDistr)
library(ggplot2)
library(bayesplot)
library(gridExtra)

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


##### specify model code #####
model_byspecies_q.file <- "model_byspecies_q.txt"
jagsscript.byspecies_q <- cat("
  model{
  # define prior: probability of q (prob of success for each sire)
  q ~ dbeta(1,1)
  
  # likelihood: stochastic relationship
  # mates ~ NB
  for(i in 1:N){
    failed.m[i] ~ dnegbin(q, nsires[i])
  }
  
  # mates ~ (negative binomial); 
  ##
  # NOTE: in R, representation of NB dist is same as Casella-Berger stat theory book
  # whereas, Ash uses the alternative NB pmf seen in the note for NB dist in Casella-Berger
  # Therefore, P(X=x|r,p) used in R is P(X=M-r|r,q) in terms of number of failures M-r
  # mean of NB is μ = n(1-p)/p and variance n(1-p)/p^2, where p is the prob of SUCCESS
  # R is calculating the number of failures M-r, and we generated the number of successes (sires, R)
  # therefore the number of mates M = failures + successes = failed.m + nsires
  
  # DERIVED QUANTITIES
  # Calculate whether multiple paternity occurred (deterministic relation)
  for(i in 1:N){
    multp[i] <- nsires[i]>1
    nmates[i] <- failed.m[i] + nsires[i]
  }
}
",file = model_byspecies_q.file)



##### read and prepare data #####
dat <- read.csv("paternity_mammals.csv")

# expand data to create multiple broods within each species (sample size * 10)
for(i in 1:nrow(dat)){
  temp.sp <- expand.grid(species = rep(dat[i,]$species, dat[i,]$nbrood))
  if(i == 1){
    temp <- temp.sp
  } else {
    temp <- rbind(temp, temp.sp)
  } 
}
mates.df <- merge(temp, dat, all.x = TRUE, all.y = FALSE)

# Remove any datapoints with NA for avgsire and avgbrood, since Poisson data gen will fail
mates.dat <- mates.df[!is.na(mates.df$avgsire),]

set.seed(36849)


##### MCMC dimensions #####
ni = 10000
nb = 1000
nc = 2 # needs to match number of initial values proposed
nt = 2
n.iter = ni + nb
  
##### parameters to monitor #####
params = c("q", "failed.m", "nsires", "multp", "nmates")

##### CREATE LIST TO HOLD ALL RESULTS FROM JAGS
jags.results <- list()

##### run the model in JAGS #####
# Initial conditions for loop
# Likelihood needs to change for each nsire value, so q calculated for each # nsire
# Start JAGS
c <- 1
lambdaS <- c(1.2,3.2) # sensible range for mammals
#starttime <- Sys.time()
for(r in lambdaS){
  
  # Bayesian estimation of q for each S and K (simulated data based loosely on real data)
  
  b <- 1
  lambdaFM <- c(1,2,3,4,5) # sensible range for mammals
  
  # Loop for mates
  for (m in lambdaFM) {
    
    a <- 1
    
    lambdaK <- c(1:10) # sensible range for mammals c(r:10)
    
    # Loop for littersize
    for (k in lambdaK) {
      
      #### generate broods and sires from Poisson for each species ####
      mates <- data.frame(cbind(rep(r, 100), rep(m, 100), rep(k, 100), matrix(NA, 100, 3)))
      names(mates) <- c("lambdaS", "lambdaFM", "lambdaK", "littersize", "nsires", "failed.m")
      for(i in 1:nrow(mates)){
        # generate broods
        avg.bs <- mates[i,]$lambdaK  # average brood size
        gen.k <- rtpois(1, avg.bs, a = 1, b = Inf)
        mates[i,]$littersize <- gen.k
  
        # generate sires
        avg.sire <- mates[i,]$lambdaS  # average number of sires
        sire.max <- mates[i,]$littersize # max # sires per litter is litter size

        gen.r <- rtpois(1, avg.sire, a = 0, b = sire.max)
        mates[i,]$nsires <- gen.r
  
        # generate failed mates
        mates[i,]$failed.m <- mates[i,]$lambdaFM
      }
    
      m.dat <- mates

      # data containing only relevant values
      jags.dat <- list(nsires = m.dat$nsires,
                       failed.m = m.dat$failed.m,
                       N = nrow(m.dat))
      
      # run JAGS (not using initial values for q)
      jmod <- jags.model(file = model_byspecies_q.file, 
                         data = jags.dat, n.chains = nc, n.adapt = 1000)
      update(jmod, n.iter = nb, by = 1, progress.bar = 'text')
      post <- coda.samples(jmod, params, n.iter = ni, thin = nt) 
      
      # save current JAGS output
      name <- paste0("lambdaS=",r," lambdaFM=",m," lambdaK=",k)
      jags.results[[name]] <- post
      
      # Save quantities from data
      MCMC_temp.k <- cbind(
        as.numeric(r),
        as.numeric(m),
        as.numeric(k),
        post.summ(post, "q")[[1]],  # mean of the est. param. q
        post.summ(post, "q")[[2]],  # sd of the est. param. q
        post.summ(post, "q")[[4]],  # 2.5% est. param. q
        post.summ(post, "q")[[5]],  # 97.5% est. param. q
        post.summ(post, "failed.m")[[1]],  # mean of failed.m
        post.summ(post, "failed.m")[[2]],  # sd of failed.m
        post.summ(post, "failed.m")[[4]],  # 2.5% failed.m
        post.summ(post, "failed.m")[[5]],  # 97.5% failed.m
        post.summ(post, "nsires")[[1]],  # mean of nsires
        post.summ(post, "nsires")[[2]],  # sd of nsires
        post.summ(post, "nsires")[[4]],  # 2.5% nsires
        post.summ(post, "nsires")[[5]],  # 97.5% nsires
        post.summ(post, "multp")[[1]],  # mean of multp
        post.summ(post, "multp")[[2]],  # sd of multp
        post.summ(post, "multp")[[4]],  # 2.5% multp
        post.summ(post, "multp")[[5]],  # 97.5% multp
        mean(m.dat$littersize),  # mean of littersize (should match lambdaK)
        var(m.dat$littersize),  # var of littersize
        post.summ(post, "nmates")[[1]],  # mean of nmates
        post.summ(post, "nmates")[[2]]  # sd of nmates
        )
      
      if(a==1) MCMC_sumtemp.k <- MCMC_temp.k else MCMC_sumtemp.k <- rbind(MCMC_sumtemp.k, MCMC_temp.k)
      
      # move to next step
      a <- a + 1
    }
    
    if(b==1) MCMC_sumtemp.FM <- MCMC_sumtemp.k else MCMC_sumtemp.FM <- rbind(MCMC_sumtemp.FM, MCMC_sumtemp.k)
    
    # move to next step
    b <- b + 1
    
  }
  
  if(c==1) MCMC_sumtemp <- MCMC_sumtemp.FM else MCMC_sumtemp <- rbind(MCMC_sumtemp, MCMC_sumtemp.FM)
  
  # move to next step
  c <- c + 1
  
}


# Species for each k
MCMC_summary <- data.frame(lambdaS = as.numeric(MCMC_sumtemp[,1]),
                           lambdaFM = as.numeric(MCMC_sumtemp[,2]),
                           lambdaK = as.numeric(MCMC_sumtemp[,3]),
                           mean.q = as.numeric(MCMC_sumtemp[,4]),
                           sd.q = as.numeric(MCMC_sumtemp[,5]),
                           q2.5 = as.numeric(MCMC_sumtemp[,6]),
                           q97.5 = as.numeric(MCMC_sumtemp[,7]),
                           mean.failedm = as.numeric(MCMC_sumtemp[,8]),
                           sd.failedm = as.numeric(MCMC_sumtemp[,9]),
                           failedm2.5 = as.numeric(MCMC_sumtemp[,10]),
                           failedm97.5 = as.numeric(MCMC_sumtemp[,11]),
                           mean.nsires = as.numeric(MCMC_sumtemp[,12]),
                           sd.nsires = as.numeric(MCMC_sumtemp[,13]),
                           nsires2.5 = as.numeric(MCMC_sumtemp[,14]),
                           nsires97.5 = as.numeric(MCMC_sumtemp[,15]),
                           mean.multp = as.numeric(MCMC_sumtemp[,16]),
                           sd.multp = as.numeric(MCMC_sumtemp[,17]),
                           multp2.5 = as.numeric(MCMC_sumtemp[,18]),
                           multp97.5 = as.numeric(MCMC_sumtemp[,19]),
                           mean.littersize = as.numeric(MCMC_sumtemp[,20]),
                           var.littersize = as.numeric(MCMC_sumtemp[,21]),
                           mean.nmates = as.numeric(MCMC_sumtemp[,22]),
                           sd.nmates = as.numeric(MCMC_sumtemp[,23]))


#### q vs k varying S plot
q_summary <- MCMC_summary
facetlabs <- c("Failed mates = 1","Failed mates = 2","Failed mates = 3","Failed mates = 4","Failed mates = 5")
q_summary$lambdaFM <- as.factor(q_summary$lambdaFM)
levels(q_summary$lambdaFM) <- facetlabs

minS1 <- min(q_summary[q_summary$lambdaS==1.2,]$mean.nsires)
maxS1 <- max(q_summary[q_summary$lambdaS==1.2,]$mean.nsires)
minS2 <- min(q_summary[q_summary$lambdaS==3.2,]$mean.nsires)
maxS2 <- max(q_summary[q_summary$lambdaS==3.2,]$mean.nsires)
q_summary$lambdaS <- as.factor(q_summary$lambdaS)
linelabs <- c(paste0(minS1, " \u2264 S \u2264 ",  maxS1), paste0(minS2, " \u2264 S \u2264 ",  maxS2))
  
qvsk_varyS <- ggplot(q_summary[q_summary$lambdaFM=="Failed mates = 3",], 
                     aes(group = as.factor(lambdaS), 
                         linetype = as.factor(lambdaS))) + 
  geom_point(aes(x=mean.littersize, y=mean.q)) + 
  geom_line(stat = "smooth", method = "loess",
            aes(x=mean.littersize, y=mean.q)) +
  geom_errorbar(aes(x=mean.littersize,
                    ymin=mean.q-(1.65*sd.q), ymax=mean.q+(1.65*sd.q)), width=.1) +
  theme_bw() + 
  scale_linetype_manual("Average # sires", values = as.factor(lambdaS), labels = linelabs) +
  ylim(c(0.2,0.5)) +
  xlab("mean littersize") + 
  ylab("q") +
  theme(legend.key.width=unit(2,"line"), # widen legend to see line type
        legend.position="bottom") 

qvsk_varyS

# ggsave("qvsk_varyS_sim.eps", plot = qvsk_varyS, device = cairo_ps,
#        scale = 1, width = 8, height = 6, units = "in",
#        dpi = 320)


#### q vs k varying FM plot ####
qvsk_varyFM <- ggplot(q_summary[q_summary$lambdaS==3.2,], 
                     aes(group = as.factor(lambdaFM), 
                         linetype = as.factor(lambdaFM))) + 
  geom_point(aes(x=mean.littersize, y=mean.q)) + 
  geom_line(stat = "smooth", method = "loess",
            aes(x=mean.littersize, y=mean.q)) +
  geom_errorbar(aes(x=mean.littersize,
                    ymin=mean.q-(1.65*sd.q), ymax=mean.q+(1.65*sd.q)), width=.1) +
  theme_bw() + 
  scale_linetype_manual("Average # failed mates", values = c(2,3,4,1,5), labels = facetlabs) +
  ylim(c(0.2,0.85)) +
  xlab("mean littersize") + 
  ylab("q") +
  theme(legend.key.width=unit(2,"line"), # widen legend to see line type
        legend.position="bottom") 

qvsk_varyFM

# ggsave("qvsk_varyFM_sim.eps", plot = qvsk_varyFM, device = cairo_ps,
#        scale = 1, width = 10, height = 6.5, units = "in",
#        dpi = 320)

