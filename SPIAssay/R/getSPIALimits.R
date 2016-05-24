######################################################################
#
# getSPIALimits: computes the limits for the SPIA test
# 
# Input: 1. N: number of valid calls
#        2. Pmm: probability of mismatch in a matching population (dafault 0.1)
#        3. nsigma: parameter that characterize the limit for mm (dafault 2)
#        4. Pmm_nonM: probability of mismach in a non matching population (dafault 0.6)
#        5. nsigma_nonM: parameter that characterize the limit for mm_nonM (dafault 3)
#        6. verbose: print verbose information on error
#
######################################################################

getSPIALimits<-function(N,Pmm = 0.1, nsigma = 2, Pmm_nonM = 0.6, nsigma_nonM = 3, verbose = FALSE){
  
  #Pmm has to be less than Pmm_nonM 
  if (Pmm >= Pmm_nonM){
    if (verbose){
      message("SPIA: error in function getSPIALimits the probability of mismatch in a matching population has to be less than the mismatch probability in a mismatching population")
    }
    return(-1)
  }
  
  ##Gaussian approximationfor the binomial distribution B(N,Pmm)
  sigmasq_mm<-(N*Pmm*(1-Pmm)); 
  media_mm<-N*Pmm;   
  
  ##Gaussian approximationfor the binomial distribution B(N,Pmm_nonM)
  sigmasq_mm_nonM<-(N*Pmm_nonM*(1-Pmm_nonM)); 
  media_mm_nonM<-N*Pmm_nonM;
  
  ##Compute limits
  k.inf <- media_mm+(nsigma*sqrt(sigmasq_mm))  
  k.sup <- max(media_mm_nonM-(nsigma_nonM*sqrt(sigmasq_mm_nonM)),0)
  
  
  if (k.inf > k.sup){
    if (verbose){
      message("SPIA: error in function getSPIALimits the inf limit is greater than the sup limit.")
      message("SPIA: verify your input parameters:")
      message("SPIA: P mismatch in a matching population:", Pmm)    
      message("SPIA: P mismatch in a mismatching population:", Pmm_nonM)
    }
    return(list(liminf=k.inf/N,limsup=k.sup/N,err=TRUE ))
  }
    
  return(list(liminf=k.inf/N,limsup=k.sup/N,err=FALSE ))
}
