library(TMB)
library(covafillr)
system("touch ssmtest.cpp")
compile("ssmtest.cpp","-g -O0",
        CXXFLAGS = cxxFlags())

dyn.load(dynlib("ssmtest"))

coord <- as.matrix(expand.grid(-5:5,-5:5))
covObs <- apply(coord,1,function(x)0.1 * sum(x^2) + rnorm(1,0,0.1))

set.seed(123)
## Simulate track
n <- 30
x <- array(0,dim=c(2,n))
for(i in 2:n)
    x[,i] <- x[,i-1] + rnorm(2,0,0.1)
## Simulate obs
obs <- array(NA,dim=c(3,n))
for(i in 1:n){
    obs[1:2,i] <- x[,i] + rnorm(2,0,0.2)
    obs[3,i] <- 0.1 * sum(x[,i]^2) + rnorm(1,0,0.05)
}
    
##plot.ts(t(obs))

dat <- list(obs = obs,
            coord = coord,
            covObs = covObs,
            p = 2,
            h = c(5,5))
param <- list(logObsSd = 0,
              logObsTSd = 10,
              logStatSd = 0,
              x = matrix(0,2,n))

obj <- MakeADFun(data = dat,
                 parameters = param,
                 random = "x",
                 DLL = "ssmtest")

obj$fn()
system.time((obj$fn()))
obj$gr()

opt <- nlminb(obj$par,obj$fn,obj$gr)
