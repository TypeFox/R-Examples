###################################################################
## svdNorm:  a function to create theta surrogates from 
## 0/1 item response matrix 
svdNorm <- function(data){  
 
  #Compute theta surrogates: thetaInit
  svdOUT <- svd(sweep(data,2,apply(data,2,mean)))

  ## needed to deal with round off error
  svdu <- round(svdOUT$u[, 1], 10) 
   
  #if sum of 1st right singular vector < 0 then 
  #multiply 1st left singular vector by -1
  if(sum(svdOUT$v[,1]) < 0) svdu<- -svdu
 
  # convert into ranks
  u1.ranks <- rank(svdu)
 
  # convert into percentiles
  Ptiles <-u1.ranks/length(u1.ranks)
  
  # to avoid infinity
  Ptiles[Ptiles==1] <- 1 - 1e-10
 
  # convert into quantiles of N(0,1)
  qnorm(Ptiles, mean=0, sd = 1)
}  

