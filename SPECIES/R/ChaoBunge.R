
ChaoBunge=function(n,t=10,conf=0.95){

  ## ===========================================================================================
  ## Purpose: This function calculates the coverage-duplication estimator by Chao and Bunge 2002
  ## input:   n --- frequency and frequency data n, matrix format, numeral data frame
  ##          t --- integer cutoff value defining less abundant species, default is 10
  ##          conf--confidence level, a numerical value<1, default .95
  ## output:  coverage-duplication estimator, standard error and  confidence interval.
  ## ===========================================================================================


  if (t!=round(t)||t<0) stop("Error: The cutoff t to define less abundant species must be non-negative integer!")
  if (ncol(n)!=2||is.numeric(n[,1])==FALSE||is.numeric(n[,2])==FALSE) {
    stop("Error: The frequency of frequencies data n must be a 2-column  matrix or data frame!")
  }
  if(is.numeric(conf)==FALSE||conf>1||conf<0) stop("Error: confidence level must be a numerical value between 0 and 1, e.g. 0.95")

  if(t>nrow(n)){
    cat("Warnings: the t that defines the abundant/rare species must be smaller than number of rows in your data matrix n!","\n")
    cat("          We use t=", nrow(n), "in this calculation!","\n\n")
    t=nrow(n)
  }
  n=as.matrix(n)
  m=max(n[,1])
  ntemp=cbind(c(1:m),rep(0,m))
  ntemp[n[,1],2]=n[,2]
  n=as.matrix(ntemp)
   colnames(n)=NULL
  rownames(n)=NULL 

  ############################
  A1=sum(n[1:t,2])
  A2=sum(n[1:t,2])-n[1,2]
  B1=sum(n[1:t,1]*n[1:t,2])
  B2=B1-n[1,2]
  C2=sum(n[1:t,1]*(n[1:t,1]-1)*n[1:t,2])
  f=n[1:t,2]
  
  ## point estimate (Equation 2, page 533 of Chao and Bunge 2002)
  if(t==nrow(n)){
    theta = 1-n[1,2]*sum(n[,1]^2*n[,2])/(sum(n[,1]*n[,2]))^2
    chao40=sum(n[(2:t),2])/theta
    chao4 =(1/theta-1)*sum(n[2:nrow(n),2])-n[1,2]+sum(n[,2])}
  if(t<nrow(n)){
    theta = 1-n[1,2]*sum(n[1:t,1]^2*n[1:t,2])/(sum(n[1:t,1]*n[1:t,2]))^2
    chao4 = sum(n[(t+1):nrow(n),2])+sum(n[(2:t),2]/theta)
    chao40=sum(n[(2:t),2])/theta
  }
  chao4
  
  ## duplication standard error(Unnumbered Eq. right below Eq. 2, page 534 of Chao and Bunge 2002)
  D=sum(n[1:t,1]*n[1:t,1]*n[1:t,2])
  partial3=numeric(t)
  partial3[1]=sum(-f[2:t]/(1-f[1]*D/B1^2)^2*(-(D+f[1])/(B1^2)+(2*f[1]*D/B1^3)))
  
  for (i in 2:(t)){
    partial3[i]=1/(1-f[1]*D/B1^2)-A2/(1-f[1]*D/B1^2)^2*(-f[1]*i^2/B1^2+2*f[1]*D*i/B1^3)
  }
    
  cova3=matrix(0,t,t)
  
  for (i in 1:t){
    for (j in 1:t){
      if(i==j){
        cova3[i,j]=f[i]*(1-f[i]/chao40)		
      }
      if(i!=j){
        cova3[i,j]=-f[i]*f[j]/chao40	
      }	
    }
  }
  
  SE=0
  for ( i in 1:t){
    for (j in 1:t){
      SE=SE+partial3[i]*partial3[j]*cova3[i,j]
    }   
  }
  SE=partial3%*%cova3%*%partial3
  SE4=sqrt(SE[1,1])
  
  ##confidence interval
  coe=qnorm(1-(1-conf)/2,0,1)
  C=exp(coe*log(1+SE4^2/(chao4-A1)^2)^0.5)
  lb=A1+(chao4-A1)/C
  ub= A1+(chao4-A1)*C
  #cat("Method: Chao and Bunge 2002 coverage-duplication method ","\n\n")
  #cat("N-estimate will be calculated using less abundant species with frequency <= ",t,"!","\n\n")
  CI0=matrix(c(lb,ub),1,2)
  colnames(CI0)=c("lb","ub")
  return(list(Nhat=round(chao4),SE=SE4,CI=round(CI0)))
}




