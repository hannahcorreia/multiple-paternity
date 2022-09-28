<!--- Readme for Correia et al. 2021 in BioEssays --->
# Multiple paternity and the number of offspring
[![Publication Status: Published](https://img.shields.io/badge/Publication%20Status-Published-success)](https://doi.org/10.1002/bies.202000247)

This is a copy of the [Harvard Dataverse data repository](https://doi.org/10.7910/DVN/3DDNQZ) for publicly available code and data to reproduce analyses in  
Correia, H. E., Abebe, A., & Dobson, F. S. (2021) Multiple paternity and the number of offspring: A model reveals two major groups of species. BioEssays. 43(4). https://doi.org/10.1002/bies.202000247  

## Code
R code for analyses presented in [Correia et al. (2021)](https://doi.org/10.1002/bies.202000247) are provided in [`bioessays2_allspecies_code.R`](bioessays2_allspecies_code.R).  
JAGS must be installed on the computer/machine before being able to run this code. Please visit https://mcmc-jags.sourceforge.io/ for more info.

## Data
Data used in the analyses and required for the R code are provided:  
- [paternity_fish.csv](paternity_fish.csv)  
- [paternity_herps.csv](paternity_herps.csv)  
- [paternity_inverts.csv](paternity_inverts.csv)  
- [paternity_mammals_Avise.csv](paternity_mammals_Avise.csv)  

Variables are as follows (further information can be found in [Correia et al. 2021](https://doi.org/10.1002/bies.202000247)):  

`nbrood` - number of litters for which multiple paternity was determined for each population  
`avgbrood` - average litter size (number of offspring per litter) of each population  
`minbood` - minimum litter size  
`maxbrood` - maximum litter size  
`pmult` - probability of multiple paternity determined by analyses of microsatellite DNA  
`avgsire` - average number of sires per litter of each population  
`minsire` - minimum number of sires  
`maxsire` - maximum number of sires  


All data are provided herein are previously published and available from the following articles:  
- Avise, J.C., and Liu, J.X. (2010). Multiple mating and its relationship to alternative modes of gestation in male-pregnant versus female-pregnant fish species. Proceedings of the National Academy of Sciences of the USA, 107, 18915–18920.  
- Avise, J.C., and Liu, J.X. (2011). Multiple mating and its relationship to brood size in pregnant fishes versus pregnant mammals and other vivparous vertebrates. Proceedings of the National Academy of Sciences of the USA, 108, 7091–7095.  
- Avise, J.C., Tatarenkov, A., and Liu, J.X. (2011). Multiple mating and clutch size in invertebrate brooders versus pregnant vertebrates. Proceedings of the National Academy of Sciences of the USA, 108, 11512–11517.  

\
\
_By using the contents on this Github repo and the published article, you agree to cite:_

Correia, H. E., Abebe, A., & Dobson, F. S. (2021) Multiple paternity and the number of offspring: A model reveals two major groups of species. BioEssays. 43(4). https://doi.org/10.1002/bies.202000247

_and_

Correia, Hannah; Dobson, F. Stephen; Abebe, Asheber, 2022, "Replication Data for: "Multiple paternity and the number of offspring: A model reveals two major groups of species"", [https://doi.org/10.7910/DVN/3DDNQZ](https://doi.org/10.7910/DVN/3DDNQZ), Harvard Dataverse, V1, UNF:6:4Mi7J6WAWWM+ujS9uEw+gA==[fileUNF] 
