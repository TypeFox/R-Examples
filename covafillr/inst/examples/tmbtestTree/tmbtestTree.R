library(TMB)
library(covafillr)

compile("tmbtestTree.cpp",CXXFLAGS = cxxFlags())

dyn.load(dynlib("tmbtestTree"))

coord <- as.matrix(expand.grid(seq(-1,1,len=101),seq(-1,1,len=101)))
ftrue <- function(x)sum(x^2)
covObs <- apply(coord,1,function(x)ftrue(x) )##+ rnorm(1,0,0.0))

dat <- list(coord = coord,
            covObs = covObs,
            p = 2,
            h = c(2,2),
            d = 100)

system.time(obj <- MakeADFun(data = dat,
                 parameters = list(x = c(0.5,0.5)),
                 DLL = "tmbtestTree"))


obj$fn(c(3.2,5.4))
obj$fn(c(0,0))
obj$fn(c(0,1))
obj$gr()

fn <- Vectorize(function(x,y)obj$fn(c(x,y)))

system.time(fn(2,6))

x <- y <- seq(-1,1,len=201)

system.time(ztrue <- outer(x,y,Vectorize(function(x,y)ftrue(c(x,y)))))

system.time(zfit <- outer(x,y,fn))

lfit0 <- loess(covObs ~ x+y,data=data.frame(covObs=covObs,x=coord[,1],y=coord[,2]))
system.time(lfit <- outer(x,y,function(x,y)predict(lfit0,newdata=data.frame(x=x,y=y))))

dev.new()
par(mfrow=c(1,3))
image(x,y,ztrue)
image(x,y,zfit)
image(x,y,lfit)

####################
## TEST MED DIM 1 ##
####################

dyn.load(dynlib("tmbtestTree"))

##coord <- as.matrix(expand.grid(seq(-1,1,len=101),seq(-1,1,len=101)))
coord <- as.matrix(seq(0,1,len=101))
ftrue <- function(x)cos(x)
grtrue <- function(x)-sin(x)
covObs <- apply(coord,1,function(x)ftrue(x) )##+ rnorm(1,0,0.0))

dat <- list(coord = coord,
            covObs = covObs,
            p = 2,
            h = c(0.5),
            d = 10)

system.time(obj <- MakeADFun(data = dat,
                 parameters = list(x = c(0.5)),
                 DLL = "tmbtestTree"))

fn <- Vectorize(function(x)obj$fn(c(x)))
gr <- Vectorize(function(x)obj$gr(c(x)))

x <- seq(0,1,len=201)
ztrue <- sapply(x,ftrue)
zfit <- sapply(x,fn)
grtrue <- sapply(x,grtrue)
grfit <- sapply(x,gr)

dev.new()
par(mfrow=c(1,2))
plot(x,ztrue)
lines(x,zfit,col='red')
plot(x,grtrue)
lines(x,grfit,col='red')
