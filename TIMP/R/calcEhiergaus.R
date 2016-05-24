"calcEhiergaus" <-
function(theta, lambda, nupower)
{
  spec <- matrix(data=0,ncol = length(theta), nrow = length(lambda))
  for(i in 1:length(theta)) {
    npare<-ifelse(length(theta[[i]])==3, 3, 4) 
    for(j in 1:(length(theta[[i]])/npare)){
      joff <- (j - 1) * npare
      a<-ifelse(npare==3,1,theta[[i]][joff + 4])
      if(theta[[i]][joff+3]!=0){
        spec[, i] <- spec[, i] + 
          a * skew(theta[[i]][joff + 1], 
                   theta[[i]][joff + 2], 
                   theta[[i]][joff + 3], l2nu(lambda), 
                   nupower)
      }
      else
        spec[, i] <- spec[, i] + a * gaus(theta[[i]][joff + 1], 
                                         theta[[i]][joff + 2], 
                                         l2nu(lambda), nupower)
      
    }
    
  }
  spec 
}
