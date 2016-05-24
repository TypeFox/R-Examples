
##################################################################################################

ChaoLee1992=function(n,t=10,method="all",conf=0.95){

  ## =================================================================================
  ## Purpose: This function calculates the lower bound estimator by chao 1984
  ## input:   n --- frequency and frequency data n, matrix format, numeral data frame
  ##          t --- integer cutoff value defining less abundant species, default is 10
  ##          method --- "ACE" or "ACE-1" or "all".
  ##          conf--confidence level, a numerical value<1, default .95
  ## output:  ACE,ACE-1 estimators, standard errors, and confidence interval.
  ## =================================================================================


  if (t!=round(t)||t<0) stop("Error: The cutoff t to define less abundant species must be non-negative integer!")
  if (ncol(n)!=2||is.numeric(n[,1])==FALSE||is.numeric(n[,2])==FALSE) {
    stop("Error: The frequency of frequencies data n must be a 2-column  matrix.")
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
  estimate=numeric(4)
  standError=numeric(4)
  A1=sum(n[1:t,2])
  A2=sum(n[1:t,2])-n[1,2]
  B1=sum(n[1:t,1]*n[1:t,2])
  B2=B1-n[1,2]
  C2=sum(n[1:t,1]*(n[1:t,1]-1)*n[1:t,2])
  f=n[1:t,2]
  
  ## Point estimator (Eq. (2.14) for Chao2 (ACE) and Eq. (2.15) for Chao3 (ACE-1) on page 211, Chao and Lee 1992)
  total = sum(n[1:t,1]*n[1:t,2])
  C = 1-n[1,2]/total
  good = sum(n[1:t,2])/C
  
  gama = max(sum(good*n[1:t,1]*(n[1:t,1]-1)*n[1:t,2])/total/(total-1)-1,0)
  gama2 =gama*(1+total*(1-C)*sum(n[1:t,1]*(n[1:t,1]-1)*n[1:t,2]/total/(total-1)/C))
  chao20=good+n[1,2]/C*gama
  chao30=good+total*(1-C)/C*gama2
  

  if(t!=nrow(n)){
    chao2 = good+n[1,2]/C*gama +sum(n[(t+1):nrow(n),2])
  }
  if(t==nrow(n)){
    chao2 = good+n[1,2]/C*gama
  }  
  
  if(t!=nrow(n)){
    chao3 = good+total*(1-C)/C*gama2+sum(n[(t+1):nrow(n),2])
  }
  
  if(t==nrow(n)){
    chao3 = good+total*(1-C)/C*gama2
  }
  
  ## standard error for chao2 (Eq. (2.17) for Chao2 (ACE) estimator on page 211, Chao and Lee 1992.
  ## Formula for Chao3 not provided, but similar

  partial=numeric(t)
  partial[1]=A2/B2+C2/B2^2*((A1*B1+f[1]*B1+f[1]*A1)/(B1-1)-f[1]*A1*B1/(B1-1)^2)
  for (i in 2:(t)){
    partial[i]=(B1+i*A2)/B2-i*A2*B1/B2^2+(f[1]*C2*B1+i*(i-1)*f[1]*A1*B1+i*f[1]*A1*C2)/(B2^2*(B1-1))-(f[1]*A1*B1*C2*(B2^2*i+2*i*(B1-1)*B2))/(B2^2*(B1-1))^2
  }
  
  cova=matrix(0,t,t)
  
  for (i in 1:t){
    for (j in 1:t){
      if(i==j){ 
        cova[i,j]=f[i]*(1-f[i]/chao20)		 
      } 
      if(i!=j){
        cova[i,j]=-f[i]*f[j]/chao20	
      }	
    }
  } 
  
  SE=0
  for ( i in 1:t){
    for (j in 1:t){
      SE=SE+partial[i]*partial[j]*cova[i,j]
    }   
  }
  SE=partial%*%cova%*%partial
  SE2=sqrt(SE[1,1])

  
  ## standard error for chao3 (Formula for standard error of Chao3 was not provided by Chao and Lee 1992, but similar to Chao2).
 
  partial2=partial
  partial2[1]=partial2[1]+(2*f[1]*A1*B1*C2^2+f[1]^2*B1*C2^2+f[1]^2*A1*C2^2)/(B2^3)/(B1-1)^2-(2*f[1]^2*A1*B1*C2^2)/(B2^3*(B1-1)^3)-(2*f[1]*B1*C2+f[1]^2*C2)/(B2^2*(B1-1))+(f[1]^2*B1*C2)/(B2^2*(B1-1)^2)
  
  for (i in 2:(t)){  partial2[i]=partial2[i]+(f[1]^2*B1*C2^2+f[1]^2*A1*i*C2^2+2*f[1]^2*A1*B1*C2*i*(i-1))/(B2^3*(B1-1)^2)-(f[1]^2*A1*B1*C2^2)*(3*B2^2*i*(B1-1)^2+2*(B1-1)*i*B2^3)/(B2^3*(B1-1)^2)^2-(f[1]^2*i*C2+f[1]^2*B1*i*(i-1))/(B2^2*(B1-1))+(f[1]^2*B1*C2*(2*B2*(B1-1)*i+B2^2*i))/(B2^2*(B1-1))^2
                   }
  cova2=matrix(0,t,t)
  
  for (i in 1:t){
    for (j in 1:t){
      if(i==j){
        cova2[i,j]=f[i]*(1-f[i]/chao30)		
      }
      if(i!=j){
        cova2[i,j]=-f[i]*f[j]/chao30	
      }	
    }
  }
  
  SE=0
  for ( i in 1:t){
    for (j in 1:t){
      SE=SE+partial2[i]*partial2[j]*cova2[i,j]
    }   
  }
  SE=partial2%*%cova2%*%partial2
  SE3=sqrt(SE[1,1])
  
  CI0=matrix(0,2,2)
  
  ##confidence interval
  coe=qnorm(1-(1-conf)/2,0,1)

  C=exp(coe*log(1+SE2^2/(chao2-A1)^2)^0.5)
  CI0[1,1]=A1+(chao2-A1)/C
  CI0[1,2]= A1+(chao2-A1)*C
  C=exp(coe*log(1+SE3^2/(chao3-A1)^2)^0.5)
  CI0[2,1]=A1+(chao3-A1)/C
  CI0[2,2]= A1+(chao3-A1)*C
  colnames(CI0)=c("lb","ub")
  rownames(CI0)=c("ACE","ACE-1")
  ACE=round(chao2)
  ACE1=round(chao3)
  CI0=round(CI0)
  
  if (toupper(method)=="ACE"){
    return(list(Nhat=ACE,SE=SE2,CI=matrix(CI0[1,],1,2)))
  } else if (toupper(method)=="ACE-1"){
    return(list(Nhat=ACE1,SE=SE3,CI=matrix(CI0[1,],1,2)))
  } else if (toupper(method)=="ALL"){
    CI=CI0
    return(list(Nhat=c(ACE,ACE1),SE=c(SE2,SE3),CI=CI0))
  }
  else{
    cat("The specified METHOD is invalid. It has to be 'ACE, ACE-1','all'","\n\n") 
  }
}



