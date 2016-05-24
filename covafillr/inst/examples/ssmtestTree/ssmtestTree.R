library(TMB)
library(covafillr)
system("touch ssmtest.cpp")
compile("ssmtestTree.cpp","-g -O0",
        CXXFLAGS = cxxFlags())

dyn.load(dynlib("ssmtestTree"))

coord <- as.matrix(expand.grid(-5:5,-5:5))
covObs <- apply(coord,1,function(x)0.1 * sum(x^2))

set.seed(123)
## Simulate track
n <- 500
x <- array(0,dim=c(2,n))
for(i in 2:n)
    x[,i] <- x[,i-1] + rnorm(2,0,0.1)
## Simulate obs
obs <- array(NA,dim=c(3,n))
for(i in 1:n){
    obs[1:2,i] <- x[,i] + rnorm(2,0,0.2)
    obs[3,i] <- 0.1 * sum(x[,i]^2) + rnorm(1,0,0.2)
}
    
##plot.ts(t(obs))

dat <- list(obs = obs,
            coord = coord,
            covObs = covObs,
            p = 2,
            h = c(2,2),
            d = 1.5)
param <- list(logObsSd = 0,
              logObsTSd = 2,
              logStatSd = 0,
              x = obs[1:2,])

obj <- MakeADFun(data = dat,
                 parameters = param,
                 random = "x",
                 DLL = "ssmtestTree")

obj$fn()
system.time((obj$fn()))
obj$gr()

opt <- nlminb(obj$par,obj$fn,obj$gr)
opt


numDeriv::grad(obj$fn,obj$env$last.par[-obj$env$random])
obj$gr(obj$env$last.par[-obj$env$random])
