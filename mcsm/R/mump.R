mump=function(Nsim=10^4,MM=10){
#multiple chains for pump failure data

data=c(5,1,5,14,3,19,1,1,4,22)
time=c(94.32,15.72,62.88,125.76,5.24,31.44,1.05,1.05,2.10,10.48)
alpha=1.8;gamma=.01;delta=1

multb=NULL;M=10
for (m in 1:MM){

  lambda=jitter(data/time)*rexp(M)
  beta=rgamma(1,gamma+10*alpha)/(delta+sum(lambda))

  for (t in 2:Nsim){
    lambda=rbind(lambda,rgamma(10,data+alpha)/(time+beta[t-1]))
    beta=c(beta,rgamma(1,gamma+10*alpha)/(delta+sum(lambda[t])))
    }

  multb=cbind(multb,beta)
  }

library(coda)
gelman.plot(mcmc.list(mcmc(multb[,1]),mcmc(multb[,2]),mcmc(multb[,3]),
mcmc(multb[,4]),mcmc(multb[,5]),mcmc(multb[,6]),mcmc(multb[,7]),
mcmc(multb[,8]),mcmc(multb[,9]),mcmc(multb[,10])),lwd=2)

autocorr.diag(mcmc(multb))
}
