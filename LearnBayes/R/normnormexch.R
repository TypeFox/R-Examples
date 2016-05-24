normnormexch=function(theta,data){
  y=data[,1]
  sigma2=data[,2]
  mu=theta[1]
  tau=exp(theta[2])
  
  logf=function(mu,tau,y,sigma2)
    dnorm(y,mu,sqrt(sigma2+tau^2),log=TRUE)
  
  sum(logf(mu,tau,y,sigma2))+log(tau)
}