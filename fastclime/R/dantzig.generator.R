
#---------------------------------------------------------------------------------#
# Package: fastclime                                                              #
# Dantzig.generator: Generates sparse linear regression model for testing dantzig #
# Authors: Haotian Pang, Di Qi, Han Liu and Robert Vanderbei                      #
# Emails: <hpang@princeton.edu>, <hanliu@princeton.edu> and <rvdb@princetonedu>   #
# Date: April 22th 2016                                                           #
# Version: 1.4.1                                                                  #
#---------------------------------------------------------------------------------#

dantzig.generator <- function(n = 50, d = 100, sparsity = 0.1, sigma0=1)
{	
  
  if(sparsity<1){      
    s<-floor(d*sparsity)
  }
  else{
    s<-sparsity
  } 
 
        
  BETA<-rep(0,d)
  pos<-rep(0,s)

  for (i in 1:s)
  {

    a<-rnorm(1, mean=0, sd=1)
    si<-2*(rbinom(1,1,0.5)-0.5)
    n1<-floor(runif(1,min=1,max=d+1))
    BETA[n1]=si*(1+a)
    pos[i]=n1
  }

  sigma<-rnorm(n, mean=0, sd=sigma0)


  X0<-matrix(rnorm(n*d, mean = 0, sd = 1), n,d)
  y<-X0%*%BETA+sigma
      
	gc()
	
	sim = list(X0 = X0 , y = y, BETA=BETA, s=s, pos=pos)
	class(sim) = "sim" 
	return(sim)
}

