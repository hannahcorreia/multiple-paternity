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

species <- c("Agile antechinus, Antechinus agilis",
             "Alpine marmot, Marmota marmota",
             "Brown bear, Ursus arctos",
             "Columbian ground squirrel, Urocitellus columbianus",
             "Eastern chipmunk, Tamias straitus",
             "Red-backed vole, Myodes rutilus",
             "Red squirrel, Tamiasciurus hudsonicus",
             "Stuart's antechinus, Antechinus stuartii",
             "Virginia opossum, Didelphis virginiana",
             "Wild boar, Sus scrofa")

b <- 1
# Start JAGS
#starttime <- Sys.time()
for(r in species){
  
  mates.p <- mates.dat[mates.dat$species==r,]
  
  if(r==species[1]){
    lambdaFM <- c(1:5) # A.antechinus
  }else if(r==species[2]){
    lambdaFM <- c(1:5) # marmot
  }else if(r==species[3]){
    lambdaFM <- c(1:5) # bear
  }else if(r==species[4]){
    lambdaFM <- c(1:5) # CGS
  }else if(r==species[5]){
    lambdaFM <- c(1:5) # chipmunk
  }else if(r==species[6]){
    lambdaFM <- c(1:5) # vole
  }else if(r==species[7]){
    lambdaFM <- c(1:5) # red squirrel
  }else if(r==species[8]){
    lambdaFM <- c(1:5) # S.antechinus
  }else if(r==species[9]){
    lambdaFM <- c(1:5) # opossum
  }else if(r==species[10]){
    lambdaFM <- c(1:5) # boar
  }else{print("Species not found for lambdaFM")}
  
  # Bayesian estimation of q for each failedM (separately)
  a <- 1
  # Loop for mates
  for (m in lambdaFM) {
    
    #### generate broods and sires from Poisson for each species ####
    mates <- cbind(mates.p, littersize = NA, nsires = NA, failed.m = NA)
    for(i in 1:nrow(mates)){
      # generate broods
      nsample <- mates[i,]$nbrood   # number of broods
      avg.bs <- mates[i,]$avgbrood  # average brood size
      max.bs <- if(!is.na(mates[i,]$maxbrood)){mates[i,]$maxbrood}else if(!is.na(mates[i,]$maxbrood_lit)){mates[i,]$maxbrood_lit}else{Inf} # max brood size possible for species from the literature
      gen.k <- rtpois(1, avg.bs, a = 1, b = max.bs)
      mates[i,]$littersize <- gen.k
      
      # generate sires
      avg.sire <- mates[i,]$avgsire
      sire.max <- mates[i,]$littersize # max # sires per litter is litter size
      
      if(r==species[1]){
        # A.antechinus p=0.978, avg.sire=2.8, maxbrood_lit=NA
        p <- 0.978
        sirepois <- rtpois(1, (5.45-(1-p))/p, a = 1, b = sire.max) 
        # check for right value of x in (x-(1-p)/p)
        switchr <- rbern(1, 1-p) # pick whether sire is 1 or >1
        gen.r <- switchr + (1-switchr)*sirepois
      }else if(r==species[2]){
        # marmot p=0.265, avg.sire=1.3, maxbrood_lit=6
        p <- 0.265
        gen.r <- rtpois(1, (0.882-(1-p))/p, a = 0, b = sire.max)
      }else if(r==species[3]){
        # brown bear p=0.145, avg.sire=1.1, maxbrood_lit=4
        p <- 0.145
        gen.r <- rtpois(1, (.883-(1-p))/p, a = 0, b = sire.max) 
      }else if(r==species[4]){
        # CGS p=0.578, avg.sire=1.8, maxbrood_lit=7
        p <- 0.578
        sirepois <- rtpois(1, (1.68-(1-p))/p, a = 1, b = sire.max) 
        switchr <- rbern(1, 1-p) # pick whether sire is 1 or >1
        gen.r <- switchr + (1-switchr)*sirepois
      }else if(r==species[5]){
        # chipmunk p=0.860, avg.sire=2.5, maxbrood_lit=5
        p <- 0.860
        sirepois <- rtpois(1, (20-(1-p))/p, a = 1, b = sire.max)
        switchr <- rbern(1, 1-p) # pick whether sire is 1 or >1
        gen.r <- switchr + (1-switchr)*sirepois
      }else if(r==species[6]){
        # red squirrel p=0.727, avg.sire=2.0, maxbrood_lit=NA
        p <- 0.727
        sirepois <- rtpois(1, (2.65-(1-p))/p, a = 1, b = sire.max) 
        switchr <- rbern(1, 1-p) # pick whether sire is 1 or >1
        gen.r <- switchr + (1-switchr)*sirepois
      }else if(r==species[7]){
        # vole p=0.232, avg.sire=1.2, maxbrood_lit=9
        p <- 0.232
        gen.r <- rtpois(1, (.856-(1-p))/p, a = 0, b = sire.max) 
      }else if(r==species[8]){
        # S. antechinus p=0.923, avg.sire=2.8, maxbrood_lit=NA
        p <- 0.923
        sirepois <- rtpois(1, (5.41-(1-p))/p, a = 1, b = sire.max) 
        switchr <- rbern(1, 1-p) # pick whether sire is 1 or >1
        gen.r <- switchr + (1-switchr)*sirepois
      }else if(r==species[9]){
        # opossum p=0.406, avg.sire=1.4, maxbrood=8
        p <- 0.406
        gen.r <- rtpois(1, (.885-(1-p))/p, a = 0, b = sire.max) 
      }else if(r==species[10]){
        # boar p=0.100, avg.sire=1.1, maxbrood=7
        p <- 0.100
        gen.r <- rtpois(1, (.901-(1-p))/p, a = 0, b = sire.max) 
      }else {print("Species not found for lambdaFM")}
      
      mates[i,]$nsires <- gen.r
      
      # generate failed mates
      mates[i,]$failed.m <- m
    }
    
    ## remove an average sire value of NA
    m.dat <- mates[!is.na(mates$avgsire), ]
    
    # data containing only relevant values
    orig.dat <- dat[dat$species==r,]
    jags.dat <- list(nsires = m.dat$nsires,
                     failed.m = m.dat$failed.m,
                     N = nrow(m.dat))
    
    # run JAGS (not using initial values for q)
    jmod <- jags.model(file = model_byspecies_q.file, 
                       data = jags.dat, n.chains = nc, n.adapt = 1000)
    update(jmod, n.iter = nb, by = 1, progress.bar = 'text')
    post <- coda.samples(jmod, params, n.iter = ni, thin = nt) 
    
    # save current JAGS output
    name <- paste0("species=",r," lambdaFM=",m)
    jags.results[[name]] <- post
    
    # Save quantities from data
    MCMC_temp.FM <- cbind(
      as.character(orig.dat$species),
      as.numeric(orig.dat$avgbrood),
      as.numeric(orig.dat$avgsire),
      as.numeric(orig.dat$pmult),
      as.numeric(m),
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
    
    if(a==1) MCMC_sumtemp.FM <- MCMC_temp.FM else MCMC_sumtemp.FM <- rbind(MCMC_sumtemp.FM, MCMC_temp.FM)
    
    # move to next step
    a <- a + 1
  }
  
  if(b==1) MCMC_sumtemp <- MCMC_sumtemp.FM else MCMC_sumtemp <- rbind(MCMC_sumtemp, MCMC_sumtemp.FM)
  
  # move to next step
  b <- b + 1
  
}


# Species for each k
MCMC_qsummary <- data.frame(species = as.character(MCMC_sumtemp[,1]),
                           avgbrood = as.numeric(MCMC_sumtemp[,2]),
                           avgsire = as.numeric(MCMC_sumtemp[,3]),
                           pmult = as.numeric(MCMC_sumtemp[,4]),
                           lambdaFM = as.numeric(MCMC_sumtemp[,5]),
                           mean.q = as.numeric(MCMC_sumtemp[,6]),
                           sd.q = as.numeric(MCMC_sumtemp[,7]),
                           q2.5 = as.numeric(MCMC_sumtemp[,8]),
                           q97.5 = as.numeric(MCMC_sumtemp[,9]),
                           mean.failedm = as.numeric(MCMC_sumtemp[,10]),
                           sd.failedm = as.numeric(MCMC_sumtemp[,11]),
                           failedm2.5 = as.numeric(MCMC_sumtemp[,12]),
                           failedm97.5 = as.numeric(MCMC_sumtemp[,13]),
                           mean.nsires = as.numeric(MCMC_sumtemp[,14]),
                           sd.nsires = as.numeric(MCMC_sumtemp[,15]),
                           nsires2.5 = as.numeric(MCMC_sumtemp[,16]),
                           nsires97.5 = as.numeric(MCMC_sumtemp[,17]),
                           mean.multp = as.numeric(MCMC_sumtemp[,18]),
                           sd.multp = as.numeric(MCMC_sumtemp[,19]),
                           multp2.5 = as.numeric(MCMC_sumtemp[,20]),
                           multp97.5 = as.numeric(MCMC_sumtemp[,21]),
                           mean.littersize = as.numeric(MCMC_sumtemp[,22]),
                           var.littersize = as.numeric(MCMC_sumtemp[,23]),
                           mean.nmates = as.numeric(MCMC_sumtemp[,24]),
                           sd.nmates = as.numeric(MCMC_sumtemp[,25]))



#### Estimate q for all species ####
## read and prepare data 
dat <- read.csv("paternity_mammals.csv")

# Expand data
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

species <- c("Agile antechinus, Antechinus agilis",
             "Alpine marmot, Marmota marmota",
             "Brown bear, Ursus arctos",
             "Columbian ground squirrel, Urocitellus columbianus",
             "Eastern chipmunk, Tamias straitus",
             "Red-backed vole, Myodes rutilus",
             "Red squirrel, Tamiasciurus hudsonicus",
             "Stuart's antechinus, Antechinus stuartii",
             "Virginia opossum, Didelphis virginiana",
             "Wild boar, Sus scrofa")

##### specify model code #####
model_q.file <- "model_q.txt"
jagsscript_q <- cat("
   model{
    q ~ dbeta(1,1)
    for(i in 1:N){
      nsires[i] ~ dbinom(q, littersize[i]) T(0,)
    }
    for(i in 1:N){
      nmates[i] <- nsires[i] / q
    }
}
",file = model_q.file)

set.seed(36849)

##### MCMC dimensions #####
ni = 10000
nb = 1000
nc = 2 # needs to match number of initial values proposed
nt = 2
n.iter = ni + nb

##### parameters to monitor #####
params = c("q", "nsires", "littersize", "nmates")

##### CREATE LIST TO HOLD ALL RESULTS FROM JAGS
jags.results <- list() 

##### run the model in JAGS #####
# Initial conditions for loop
b <- 1
# Start JAGS
#starttime <- Sys.time()
for(r in species){
  
  mates.p <- mates.dat[mates.dat$species==r,]
  mates.p$nsires <- NA
  mates.p$littersize <- NA
  
  for(i in 1:nrow(mates.p)){
    if(r %in% species[4]){
      nbrood <- mates.p[i,]$nbrood
      avg.bs <- mates.p[i,]$avgbrood
      avgsire <- mates.p[i,]$avgsire
      max.bs <- if(!is.na(mates.p[i,]$maxbrood)){
        mates.p[i,]$maxbrood
      }else if(!is.na(mates.p[i,]$maxbrood_lit)){
        mates.p[i,]$maxbrood_lit}else{Inf}
    }else{
      nbrood <- mates.p[i,]$nbrood
      avg.bs <- mates.p[i,]$avgbrood
      avgsire <- mates.p[i,]$avgsire
      max.bs <- if(!is.na(mates.p[i,]$maxbrood)){
        mates.p[i,]$maxbrood
      }else if(!is.na(mates.p[i,]$maxbrood_lit)){
        mates.p[i,]$maxbrood_lit}else{Inf}
    }
    
    # create simulated data based on data avg k and avg s
    # generate k first
    if(r==species[1]){
      gen.k <- rtpois(1, 9.52, a = 1, b = max.bs)
    }else if(r==species[2]){
      gen.k <- rtpois(1, 3.24, a = 1, b = max.bs)
    }else if(r==species[3]){
      gen.k <- rtpois(1, 0.01, a = 1, b = max.bs)
    }else if(r==species[4]){
      ifelse(nbrood==147, 
             gen.k <- rtpois(1, 2.19, a = 1, b = max.bs),
             gen.k <- rtpois(1, 1.818, a = 1, b = max.bs))
    }else if(r==species[5]){
      gen.k <- rtpois(1, 2.8, a = 1, b = max.bs)
    }else if(r==species[6]){
      gen.k <- rtpois(1, 2.5, a = 1, b = max.bs)
    }else if(r==species[7]){
      gen.k <- rtpois(1, 4.05, a = 1, b = max.bs)
    }else if(r==species[8]){
      gen.k <- rtpois(1, 4.32, a = 1, b = max.bs)
    }else if(r==species[9]){
      gen.k <- rtpois(1, 9.265, a = 1, b = max.bs)
    }else{
      gen.k <- rtpois(1, 4.85, a = 1, b = max.bs)
    }
    
    # set sire.max as the litter size
    mates.p[i,]$littersize <- gen.k
    sire.max <- mates.p[i,]$littersize
    
    # generate r sires next
    
    if(r==species[1]){
      # A.antechinus p=0.978, avg.sire=2.8, maxbrood_lit=NA
      p <- 0.978
      sirepois <- rtpois(1, (5.45-(1-p))/p, a = 1, b = sire.max) 
      switchr <- rbern(1, 1-p) # pick whether sire is 1 or >1
      gen.r <- switchr + (1-switchr)*sirepois
    }else if(r==species[2]){
      # marmot p=0.265, avg.sire=1.3, maxbrood_lit=6
      p <- 0.265
      gen.r <- rtpois(1, (0.882-(1-p))/p, a = 0, b = sire.max)
    }else if(r==species[3]){
      # brown bear p=0.145, avg.sire=1.1, maxbrood_lit=4
      p <- 0.145
      gen.r <- rtpois(1, (.883-(1-p))/p, a = 0, b = sire.max)
    }else if(r==species[4]){
      if(nbrood==147){ 
        # CGS p=0.578, avg.sire=1.8, maxbrood_lit=7
        p <- 0.578
        sirepois <- rtpois(1, (1.68-(1-p))/p, a = 1, b = sire.max) 
        switchr <- rbern(1, 1-p) # pick whether sire is 1 or >1
        gen.r <- switchr + (1-switchr)*sirepois
      }else{
        # CGS p=0.724, avg.sire=1.9, maxbrood_lit=7
        p <- 0.724
        sirepois <- rtpois(1, (1.9-(1-p))/p, a = 1, b = sire.max) 
        switchr <- rbern(1, 1-p) # pick whether sire is 1 or >1
        gen.r <- switchr + (1-switchr)*sirepois
      }
    }else if(r==species[5]){
      # chipmunk p=0.860, avg.sire=2.5, maxbrood_lit=5
      p <- 0.860
      sirepois <- rtpois(1, (20-(1-p))/p, a = 1, b = sire.max)
      switchr <- rbern(1, 1-p) # pick whether sire is 1 or >1
      gen.r <- switchr + (1-switchr)*sirepois
    }else if(r==species[6]){
      # red squirrel p=0.727, avg.sire=2.0, maxbrood_lit=NA
      p <- 0.727
      sirepois <- rtpois(1, (2.65-(1-p))/p, a = 1, b = sire.max) 
      switchr <- rbern(1, 1-p) # pick whether sire is 1 or >1
      gen.r <- switchr + (1-switchr)*sirepois
    }else if(r==species[7]){
      # vole p=0.232, avg.sire=1.2, maxbrood_lit=9
      p <- 0.232
      gen.r <- rtpois(1, (.856-(1-p))/p, a = 0, b = sire.max) 
    }else if(r==species[8]){
      # S. antechinus p=0.923, avg.sire=2.8, maxbrood_lit=NA
      p <- 0.923
      sirepois <- rtpois(1, (5.41-(1-p))/p, a = 1, b = sire.max) 
      switchr <- rbern(1, 1-p) # pick whether sire is 1 or >1
      gen.r <- switchr + (1-switchr)*sirepois
    }else if(r==species[9]){
      # opossum p=0.406, avg.sire=1.4, maxbrood=8
      p <- 0.406
      gen.r <- rtpois(1, (.885-(1-p))/p, a = 0, b = sire.max) 
    }else if(r==species[10]){
      # boar p=0.100, avg.sire=1.1, maxbrood=7
      p <- 0.100
      gen.r <- rtpois(1, (.901-(1-p))/p, a = 0, b = sire.max) 
    }else {print("Species not found for sires")}
    
    mates.p[i,]$nsires <- gen.r
  }
  
  # data containing only relevant values
  m.dat <- mates.p
  orig.dat <- dat[dat$species==r,]
  jags.dat <- list(nsires = m.dat$nsires, 
                   littersize = m.dat$littersize, 
                   N = nrow(m.dat))
  
  # run JAGS (not using initial values for q)
  jmod <- jags.model(file = model_q.file, data = jags.dat, n.chains = nc, n.adapt = 1000)  
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
    post.summ(post, "littersize")[[1]],  # mean of est.M
    post.summ(post, "littersize")[[2]],  # sd of est.M
    post.summ(post, "littersize")[[4]],  # 2.5% est.M
    post.summ(post, "littersize")[[5]]  # 97.5% est.M
  )
  
  if(b==1) MCMC_sumtemp <- MCMC_temp else MCMC_sumtemp <- rbind(MCMC_sumtemp, MCMC_temp)
  
  # move to next step
  b <- b + 1
}

MCMC_estq_estm <- data.frame(species = as.character(MCMC_sumtemp[,1]),
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
                             mean.littersize = as.numeric(MCMC_sumtemp[,17]),
                             sd.littersize = as.numeric(MCMC_sumtemp[,18]),
                             littersize2.5 = as.numeric(MCMC_sumtemp[,19]),
                             littersize97.5 = as.numeric(MCMC_sumtemp[,20]))



## Need resid p from previous data
qnull_resid <- read.csv("MCMC_resids.csv")

names(qnull_resid)[5:8] <- c("mean.qnull", "sd.qnull", "qnull2.5", "qnull97.5")

q_dat <- merge(MCMC_estq_estm, qnull_resid[,c(1,3,36)], by = c("species", "avgsire"), all.x = TRUE)

# All species q in one figure (blue = p-pB is positive species; red = p-pB is negative species)
ppb2 <- q_dat[q_dat$species %in% species,]
### Check pB-p
ppb2[, c("species", "resid_pmult")]

## clean LaTeX table of pB-p, mean.q (95% CI q), and mean.nmates (95% CI)
library(xtable)
xtable(ppb2[,c("species", "avgbrood", "avgsire", "pmult", "resid_pmult", "mean.q", "q2.5", "q97.5", 
               "mean.nmates", "nmates2.5", "nmates97.5")],
       digits=c(0,3,3,3,3,3,3,3,3,3,3,3))


#### 7 representative species ####
## Need null q from previous data
qnulldat <- MCMC_estq_estm[,c(1,3,5:8,13:16)]

# Rename variables not to confuse with q line data in MCMC_qsummary
names(qnulldat)[3:10] <- c("q_emp", "qsd_emp", "q2.5_emp", "q97.5_emp", "nmates_emp", "nmatessd_emp", "nmates2.5_emp", "nmates97.5_emp")

q_summary2 <- merge(MCMC_qsummary, qnulldat, by = c("species", "avgsire"), all.x = TRUE)

species7 <- c("Agile antechinus, Antechinus agilis",
              "Columbian ground squirrel, Urocitellus columbianus",
              "Eastern chipmunk, Tamias straitus",
              "Red squirrel, Tamiasciurus hudsonicus",
              "Stuart's antechinus, Antechinus stuartii",
              "Virginia opossum, Didelphis virginiana",
              "Wild boar, Sus scrofa")
q_summary7 <- q_summary2[q_summary2$species %in% species7,]
q_summary7$species <- droplevels(q_summary7$species)

q_summary7$q_zero <- NA
for(r in species7){
  mod1 <- mgcv::gam(mean.q ~ s(mean.nmates, k = 5), data = q_summary7[q_summary7$species==r,])
  x0 <- data.frame(mean.nmates = q_summary7[q_summary7$species==r,]$nmates_emp)
  y_zero <- predict(mod1, x0)
  q_summary7[q_summary7$species==r,]$q_zero <- y_zero
}

## blue = p-pB is positive species; red = p-pB is negative species
species7cols <- c("red", "black", "red", "blue", "red", "black", "black")
species7lines <- c(2,6,1,1,3,5,4)
myspecies7 <- c("Agile antechinus", "Columbian ground squirrel", 
                "Eastern chipmunk", "Red squirrel", "Stuart's antechinus", 
                "Virginia opossum", "Wild boar")

qvsm_7 <- ggplot(q_summary7, aes(group = as.factor(species), 
                                 color = as.factor(species),
                                 linetype = as.factor(species))) + 
  # estimated q
  geom_point(aes(x=mean.nmates, y=mean.q), shape = 16) + 
  geom_line(stat = "smooth", method = "gam", formula = y ~ s(x,k=5),
            aes(x=mean.nmates, y=mean.q)) +
  geom_errorbar(aes(x=mean.nmates,
                    ymin=mean.q-(1.98*sd.q), ymax=mean.q+(1.98*sd.q)), width=.1) +
  # empirical q from Bayesian analysis
  geom_point(aes(x=nmates_emp, y=q_emp, shape = "q from S~tB(K,q)", size = 2), 
             alpha = 0.1) +
  #geom_point(aes(x=nmates_emp, y=q_emp), alpha = 0.05, shape = 8, size = 3) +
  geom_point(aes(x=nmates_emp, y=q_zero, shape = "q from M~NB(S,q)", size = 2),
             alpha = 0.1) +
  geom_segment(aes(x=nmates_emp, y=q_zero, xend=-Inf, yend=q_zero), 
               alpha = 0.05) +
  geom_segment(aes(x=nmates_emp, y=q_zero, xend=nmates_emp, yend=-Inf), 
               alpha = 0.05) +
  theme_classic() + 
  scale_color_manual("Species", values = species7cols, labels = myspecies7) + 
  scale_linetype_manual("Species", values = species7lines, labels = myspecies7) +
  scale_shape_manual("Estimated q values", values = c(8,10), 
                     labels = c("q from M~NB(S,q)","q from S~tB(K,q)")) +
  xlab("mean # mates") + 
  ylab("q") +
  guides(shape=guide_legend(override.aes=list(size=5))) +
  guides(size=FALSE) +
  theme(legend.key.width=unit(3,"line")) # widen legend to see line type

qvsm_7

# ggsave("qvsm_empq7.eps", plot = qvsm_7, device = cairo_ps, 
#        scale = 1, width = 11, height = 7, units = "in",
#        dpi = 320)

