dBi.Gaussian=function(X,alpha=0,sigma1,sigma2)
{
  result=0
  result[which(X<alpha)]<-dnorm(X[which(X<alpha)],alpha,sigma1)*sigma1/(sigma1+sigma2)
  result[which(X>=alpha)]<-dnorm(X[which(X>=alpha)],alpha,sigma2)*sigma2/(sigma1+sigma2)
  return(result) 
}

Bi.Gaussian.Sim=function(mode,sigma1,sigma2,n)
{
  weight = c(sigma1, sigma2)
  mu = c(0,0)
  sigma = c(sigma1,sigma2)
  z = sample(c(1,2),size=n, prob=weight,replace=TRUE)
  
  r = rnorm(n,mean=mu[z],sd=sigma[z])
  if((z==1)&(r>=0)){r=-r
  }else if((z==2)&(r<=0)){r=-r}
  return(r)
}

