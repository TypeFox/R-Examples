#Function: ghap.blockstats
#License: GPLv3 or later
#Modification date: 22 Jan 2016
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Calculate block summary statistics

ghap.blockstats <- function(hapstats, ncores = 1){
  
  #Get unique blocks
  blocks <- unique(hapstats[,c("BLOCK","CHR","BP1","BP2")])
  
  #Set of internal functions
  my.fun <- function(block){
    temp<-hapstats$FREQ[hapstats$BLOCK==block]
    exp.het <- 1-(sum((temp)^2))
    return(c(exp.het,length(temp)))
  }
    
  #Calculation of expected heterozygosity
  temp <- unlist(mclapply(FUN = my.fun, blocks$BLOCK, mc.cores = ncores))
  blocks$EXP.H <- temp[1:length(temp) %% 2 == 1]
  blocks$N.ALLELES <- temp[1:length(temp) %% 2 == 0]
  
  #Return output
  return(blocks)
  
}

