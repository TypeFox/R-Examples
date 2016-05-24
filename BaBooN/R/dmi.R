# Data monotonicity index for missing values ###
# Version:                                 0.1
# Date:                             2015-03-01
# Author:                           F.M., T.S.

dmi <- function(Data){
  
  #Retrieving some informations
  dimD <- dim(Data)
  nmis <- sapply(Data,function(x) sum(is.na(x)))  
  obs_mono <- dimD[1] - nmis
  
  #Calculating sum of NAs that aren't monotone
  stillNA <- integer(dimD[2])
  
  for(i in 1:dimD[2]){
    stillNA[i] <- sum(is.na(Data[1:obs_mono[i],i]))
  }
  
  #Returning the monotonicity index
  return(1 - sum(stillNA)/sum(nmis))
}
