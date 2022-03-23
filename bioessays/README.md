<!--- Readme for Abebe et al. 2019 in BioEssays --->
# Multiple paternity and the number of offspring
[![Publication Status: Published](https://img.shields.io/badge/Publication%20Status-Published-success)](https://doi.org/10.1002/bies.201900016)

R code for analyses presented in [Abebe et al. (2019)](https://doi.org/10.1002/bies.201900016) are provided in [`bioessays_data_code.R`](bioessays_data_code.R) and [`bioessays_simulation_code.R`](bioessays_simulation_code.R).
JAGS must be installed on the computer/machine before being able to run this code. Please visit https://mcmc-jags.sourceforge.io/ for more info.

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
By using the contents on this Github repo and the published article, you agree to cite:

Correia, H. E., Abebe, A., & Dobson, F. S. (2021) Multiple paternity and the number of offspring: A model reveals two major groups of species. BioEssays. 43(4). https://doi.org/10.1002/bies.201900016
