"fillK" <- function (theta, kinscal, kmat, fixedkmat, kinscalspecial,
                     kinscalspecialspec, nocolsums = FALSE) 
  #fillK
{
  if(!fixedkmat){ 
    if(length(kinscalspecialspec) > 0) 
      special <- TRUE
    else special <- FALSE
    ##kmat is 3d array with elements kmat[i,j,]
    ##<- c(thetaindex, theta scal index)
    dimk <- nrow(kmat)  ## kmat is square 
    reskmat<-matrix(0, nrow=dimk, ncol=dimk)
    for(i in 1:dimk) {
      for(j in 1:dimk) {
        if( ! (special || nocolsums ))
          if(i==j)
            pl<- -1
        else	pl<- 1
        else pl <- 1
        
        if(!special) {## k is (possibly) linearly dep. on scaling
          if(kmat[i,j,1] != 0 && kmat[i,j,2] != 0)
            ## use scaling and theta
            reskmat[i,j]<-pl*theta[kmat[i,j,1]]*
            kinscal[kmat[i,j,2]]
          if(kmat[i,j,1] != 0 && kmat[i,j,2] == 0)	
            ## use theta 
            reskmat[i,j]<-pl*theta[kmat[i,j,1]]
          if(kmat[i,j,1] == 0 && kmat[i,j,2] != 0)	
            ## use scaling 
            reskmat[i,j]<-pl*kinscal[kmat[i,j,2]]
        }
        else{ ## have a function to describe k 
          if(kmat[i,j,2] != 0)
            reskmat[i,j] <- applyKinScalSpecialFun(ind=kmat[i,j,2],
                                                   theta = theta,                                    
                                                   kinscalspecial,kinscalspecialspec) 
          
        }
        
      }
    }
    if( ! (special || nocolsums) ) {
      for(i in 1:dimk) {
        for(j in 1:dimk) {
          if(i!=j)
            reskmat[j,j]<-reskmat[j,j]-reskmat[i,j]
        }
      }
    }
  }
  else { 
    diag(kmat) <- -  colSums(kmat)
    reskmat <- kmat
  } 
  reskmat
}

