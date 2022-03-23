<!--- Readme for Dobson et al. 2018 in Proceedings: Biological Sciences --->
[![Publication Status: Published](https://img.shields.io/badge/Publication%20Status-Published-success)](https://doi.org/10.1098/rspb.2018.2042)

R code for analyses presented in [Dobson et al. (2018)](https://doi.org/10.1098/rspb.2018.2042) are provided in [`PRSLB_code.R`](PRSLB_code.R).
JAGS must be installed on the computer/machine before being able to run this code. Please visit https://mcmc-jags.sourceforge.io/ for more info.

Data used in the analyses and required for the R code are given in [paternity_mammals.csv](paternity_mammals.csv).
Variables are as follows (further information can be found in [Dobson et al. 2018](https://doi.org/10.1098/rspb.2018.2042)):  

`nbrood` - number of litters for which multiple paternity was determined for each population  
`avgbrood` - average litter size (number of offspring per litter) of each population  
`minbood` - minimum litter size  
`maxbrood` - maximum litter size  
`pmult` - probability of multiple paternity determined by analyses of microsatellite DNA  
`avgsire` - average number of sires per litter of each population  
`minsire` - minimum number of sires  
`maxsire` - maximum number of sires  

Sources for the original data collected from manuscripts with full references are given in [PRSLB_DataAppendix.docx](PRSLB_DataAppendix.docx)


By using the contents on this Github repo and the published article, you agree to cite:  
Dobson, F. S., Abebe, A., Correia, H. E., Kasumo, C., & Zinner, B. (2018) Multiple paternity and number of offspring in mammals. Proceedings of the Royal Society B: Biological Sciences. 285(1891).[https://doi.org/10.1098/rspb.2018.2042](https://doi.org/10.1098/rspb.2018.2042)
