library(spatsurv)
library(sp)
library(spatstat)
library(survival)

par(mfrow=c(1,1))

set.seed(10)

n <- 100
DIST <- exponentialHaz()


OMEGA <- 1



# Generate spatially correlated survival data ... 
dat <- simsurv(X=cbind( age=runif(n,5,50),sex=rbinom(n,1,0.5),cancer=rbinom(n,1,0.2)),
                        dist=DIST,
                        omega=OMEGA,
                        mcmc.control=mcmcpars(nits=100,burn=10,thin=10))  

coords <- dat$coords
SIGMA <- dat$cov.parameters[1]
PHI <- dat$cov.parameters[2]                                          

par(mfrow=c(2,2))                                    
plot(coords,col=grey(1-dat$survtimes/max(dat$survtimes)),pch=19)

X <- as.data.frame(dat$X) # covariates

survtimes <- dat$survtimes
censtimes <- runif(n,min(survtimes),max(survtimes))                                    
survdat <- gencens(survtimes,censtimes)  


# priors
betaprior <- betapriorGauss(mean=0,sd=10)
omegaprior <- omegapriorGauss(mean=0,sd=10)
etaprior <- etapriorGauss(mean=log(c(SIGMA,PHI)),sd=c(0.3,0.3))
priors <- mcmcPriors(   betaprior=betaprior,
                        omegaprior=omegaprior,
                        etaprior=etaprior,
                        call=indepGaussianprior,
                        derivative=derivindepGaussianprior)

# create SpatialPointsDataFrame containing the covariate data and coordinates of the survival data 
spatdat <- SpatialPointsDataFrame(coords,data=as.data.frame(X))
spatdat$ss <- survdat

if(TRUE){
    ss <- survspat( formula=ss~age+sex+cancer,
                    data=spatdat,
                    dist=DIST,
                    cov.model=covmodel(model="exponential",pars=NULL),
                    mcmc.control=mcmcpars(nits=100,burn=10,thin=9),
                    priors=priors)
}

par(mfrow=c(3,3))
plot(ss$etasamp[,1],type="s", main="eta 1")
plot(ss$etasamp[,2],type="s", main="eta 2")
plot(ss$omegasamp[,1],type="s", main="omega")
plot(ss$betasamp[,1],type="s", main="beta 1")
plot(ss$betasamp[,2],type="s", main="beta 2")
plot(ss$Ysamp[,1],type="s", main="Y 1")
plot(ss$Ysamp[,floor(n/2)],type="s", main="Y floor(n/2)")
plot(ss$tarrec,type="s", main="log posterior")
plot(dat$Y,colMeans(ss$Ysamp),xlab="TRUE Y",ylab="Estimated Y")
abline(0,1)

#
#save(list=ls(),file="simout1.RData")
#
