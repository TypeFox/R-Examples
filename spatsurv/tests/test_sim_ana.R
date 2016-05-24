library(spatsurv)
library(survival)

set.seed(10)

dat <- simsurv(mcmc.control=mcmcpars(nits=100,burn=10,thin=9),dist=exponentialHaz())

X <- as.data.frame(dat$X) # covariates

survtimes <- dat$survtimes
n <- length(survtimes)
censtimes <- runif(n,min(survtimes),max(survtimes))                                    
survdat <- gencens(survtimes,censtimes)                                    

plot(survfit(survdat~1))
                
dat1 <- simsurv(mcmc.control=mcmcpars(nits=100,burn=10,thin=9),dist=weibullHaz(),omega=c(1,0.5))

dat1 <- simsurv(mcmc.control=mcmcpars(nits=100,burn=10,thin=9),dist=gompertzHaz(),omega=c(1,0.5))

dat1 <- simsurv(mcmc.control=mcmcpars(nits=100,burn=10,thin=9),dist=makehamHaz(),omega=c(1,0.5,1))

dat1 <- simsurv(mcmc.control=mcmcpars(nits=100,burn=10,thin=9),dist=tpowHaz(c(3,2)),omega=c(0.5,1.1))
                
