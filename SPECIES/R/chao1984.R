
chao1984=function(n,conf=0.95){

  ## ==========================================================================================
  ## Purpose: This function calculates the lower bound estimator by chao 1984
  ## input:   n--- frequency and frequency data n
  ##          conf--confidence level, a numerical value<1, default .95
  ## output:  lower bound estimator and standard error  and confidence interval.
  ## ==========================================================================================

  
  if (ncol(n)!=2||is.numeric(n[,1])==FALSE||is.numeric(n[,2])==FALSE) {
    stop("Error: The frequency of frequencies data n must be a 2-column  matrix.")
  }
  if(is.numeric(conf)==FALSE||conf>1||conf<0) stop("Error: confidence level must be a numerical value between 0 and 1, e.g. 0.95")

  n=as.matrix(n)
  m=max(n[,1])
  ntemp=cbind(c(1:m),rep(0,m))
  ntemp[n[,1],2]=n[,2]
  n=as.matrix(ntemp)
  colnames(n)=NULL
  rownames(n)=NULL 
    
  ##Chao1984--lower bound estimator (Eq. (6) on page 267, Chao 1984)
  chao1=sum(n[,2])+n[1,2]^2/2/n[2,2]
  chao1
    
  ##lower band estimator standard error (Eq.  right below Eq. (10) on page 786, Chao 1987) 
  coef=qnorm(1-(1-conf)/2,0,1)
  SE1=sqrt(n[2,2]*(0.25*(n[1,2]/n[2,2])^4+(n[1,2]/n[2,2])^3+0.5*(n[1,2]/n[2,2])^2))
  A1=sum(n[,2])
  C=exp(coef*log(1+SE1^2/(chao1-A1)^2)^0.5)

  lb=A1+(chao1-A1)/C
  ub=A1+(chao1-A1)*C
  CI0=matrix(c(lb,ub),1,2)
  colnames(CI0)=c("lb","ub")
  #cat("Method: Chao 1984 lower bound estimator ","\n\n")
  #cat("N-estimate will be calculated using less abundant species with frequency <= ",t,"!","\n\n")
  
  return(list(Nhat=round(chao1),SE=SE1,CI=round(CI0)))
}



