<!--- Readme for Abebe et al. 2019 in BioEssays --->
# Estimating a key parameter of mammalian mating systems
[![Publication Status: Published](https://img.shields.io/badge/Publication%20Status-Published-success)](https://doi.org/10.1002/bies.201900016)

This is a copy of the [Harvard Dataverse data repository](https://doi.org/10.7910/DVN/5BHZVE) for publicly available code and data to reproduce analyses in  
Abebe, A., Correia, H. E., & Dobson, F. S. (2019) Estimating a key parameter of mammalian mating systems: the chance of siring success for a mated male. BioEssays. 41(12). https://doi.org/10.1002/bies.201900016

## Code
R code for analyses presented in [Abebe et al. (2019)](https://doi.org/10.1002/bies.201900016) are provided in [`bioessays_data_code.R`](bioessays_data_code.R) and [`bioessays_simulation_code.R`](bioessays_simulation_code.R).
JAGS must be installed on the computer/machine before being able to run this code. Please visit https://mcmc-jags.sourceforge.io/ for more info.

## Data
Data used in the analyses and required for the R code are given in [paternity_mammals.csv](paternity_mammals.csv)
Variables are as follows (further information can be found in [Abebe et al. 2019](https://doi.org/10.1002/bies.201900016)):  

`nbrood` - number of litters for which multiple paternity was determined for each population  
`avgbrood` - average litter size (number of offspring per litter) of each population  
`minbood` - minimum litter size  
`maxbrood` - maximum litter size  
`pmult` - probability of multiple paternity determined by analyses of microsatellite DNA  
`avgsire` - average number of sires per litter of each population  
`minsire` - minimum number of sires  
`maxsire` - maximum number of sires  
`maxbrood_lit` - maximum brood size of species found in broader literature

R code [`bioessays_data_code.R`](bioessays_data_code.R) also requires results from [Dobson et al. (2018)](https://doi.org/10.1098/rspb.2018.2042), which are provided in [MCMC_resids.csv](MCMC_resids.csv). Code to replicate these results are in [multiple_paternity/ProcB](multiple-paternity/ProcB).

Sources for the original data collected from manuscripts with full references are given in [Dobson et al. (2018)](https://doi.org/10.1098/rspb.2018.2042)
\
\
\
_By using the contents on this Github repo and the published article, you agree to cite:_

Abebe, A., Correia, H. E., & Dobson, F. S. (2019) Estimating a key parameter of mammalian mating systems: the chance of siring success for a mated male. BioEssays. 41(12). https://doi.org/10.1002/bies.201900016

_and_

Correia, Hannah; Abebe, Asheber; Dobson, F. Stephen, 2022, "Replication Data for: "Estimating a key parameter of mammalian mating systems: the chance of siring success for a mated male"", [https://doi.org/10.7910/DVN/5BHZVE](https://doi.org/10.7910/DVN/5BHZVE), Harvard Dataverse, V1, UNF:6:qBdH9P9Mi3pn6ub4g4z+Iw== [fileUNF] 
