pkgname <- "acopula"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('acopula')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("copula")
### * copula

flush(stderr()); flush(stdout())

### Name: copula
### Title: Generic copulae definitions
### Aliases: copula copFGM copGumbel copNormal copPlackett copProduct
### Keywords: copula

### ** Examples

## the following gives the same definition list
copFGM()
copula("FGM")

## any list item can be modified upon function call
copPlackett(parameters=2.2,upper=10)



cleanEx()
nameEx("depfun")
### * depfun

flush(stderr()); flush(stdout())

### Name: depfun
### Title: Dependence function of Extreme-Value copula
### Aliases: depfun dep1 depGalambos depGumbel depHuslerReiss depMax
###   depTawn depCC depGCC ldepPartition3D
### Keywords: Pickands' dependence function Extreme-Value copula EV copula

### ** Examples

## the following gives the same definition list
depGumbel()
depfun("Gumbel")

## any list item can be modified upon function call
depGumbel(parameters=2.2,upper=10)

## general convex combination of 5 basic depfuns that arise from 
## partitioning method for 3 dimensions; it results in 
## (3x5)-parametric Pickand's dependence function definition list
depGCC(depfun=ldepPartition3D(), dim = 3)



cleanEx()
nameEx("generator")
### * generator

flush(stderr()); flush(stdout())

### Name: generator
### Title: Generator of Archimedean copula
### Aliases: generator genAMH genClayton genFrank genGumbel genJoe genLog
### Keywords: generator Archimedean copula

### ** Examples

## the following gives the same definition list
genGumbel()
generator("Gumbel")

## any list item can be modified upon function call
genGumbel(parameters=2.2,upper=10)



cleanEx()
nameEx("nderive")
### * nderive

flush(stderr()); flush(stdout())

### Name: nderive
### Title: Numerical derivative
### Aliases: nderive
### Keywords: derivative linear approximation

### ** Examples

##density of a bivariate Gumbel copula evaluated in point c(0.5,0.6)
nderive(fun = function(x) pCopula(x,genGumbel(),gpar=3.5), point = c(0.5,0.6),
        order = c(1,1))



cleanEx()
nameEx("nintegrate")
### * nintegrate

flush(stderr()); flush(stdout())

### Name: nintegrate
### Title: Numerical integration
### Aliases: nintegrate
### Keywords: integral trapezoid

### ** Examples

##cumulative distribution function of a bivariate normal copula 
##evaluated at point c(0.5,0.6); compare pCopula(c(0.5,0.6),cop=copNormal(),par=0.5)
nintegrate(function(x) dCopula(x,cop=copNormal(),par=0.5), 
  lower=0.001, upper=c(0.5,0.6), subdivisions=20) 



cleanEx()
nameEx("vpartition")
### * vpartition

flush(stderr()); flush(stdout())

### Name: vpartition
### Title: Vector partitioning
### Aliases: vpartition
### Keywords: split vector

### ** Examples

vpartition(1:10,c(4,5,2))



cleanEx()
nameEx("xCopula")
### * xCopula

flush(stderr()); flush(stdout())

### Name: xCopula
### Title: Archimax and generic copula distribution functions
### Aliases: cCopula dCopula eCopula gCopula pCopula qCopula rCopula
###   isCopula eCopulaArchimax eCopulaGeneric gCopulaEmpirical
###   pCopulaEmpirical rCopulaArchimax2D print.eCopulaArchimax
###   print.eCopulaGeneric print.gCopula print.isCopula
### Keywords: cumulative distribution function probability density function
###   conditional probability quantile sampling empirical copula
###   goodness-of-fit test d-increasing maximum likelihood inverse of
###   correlation coefficient

### ** Examples

## assign generator definition list with specific parameter
ge <- genGumbel(parameters=4)

## probability P(U<0.3,V<0.5)
pCopula(c(0.3,0.5),ge)  #0.2906142
## quantile q for which P(U<q,V<0.5)=0.2906142
pCopula(c(0.2906142,0.5),ge,quantile=1)  #0.3000175
pCopula(c(NA,0.5),ge,quantile=1,probability=0.2906142)
qCopula(c(0.5),quantile=1,probability=0.2906142,generator=ge)

## conditional probability P(U<0.3|V=0.5)
cCopula(c(0.3,0.5),ge,conditional.on=2)  #0.1025705
## quantile q for which conditional probability P(U<q|V=0.5)=0.1025705
cCopula(c(0.1025705,0.5),conditional.on=2,generator=ge,quantile=1)  #0.2999861
cCopula(c(NA,0.5),conditional.on=2,generator=ge,quantile=1,probability=0.1025705)
qCopula(c(0.5),quantile=1,probability=0.1025705,conditional.on=2,generator=ge)

## copula density
dCopula(c(0.3,0.5),ge) #1.083797
local({
x <- y <- seq(0,1,length.out=20)
persp(x,y,matrix(dCopula(expand.grid(x,y),ge),nrow=length(x)),r=2,zlab="density")
})

## simulate random vector
rge <- rCopula(100,dim=2,ge) 
plot(rge)
# Observe that using rCopula(100,dim=2,cop=copGumbel(parameters=4)) 
# would take much more time to sample, since numerical derivative needs to be employed. 

## --- fit copula to data set
# maximum likelihood (using density)
eCopula(rge,ge,technique="ML")
# some methods has no support for parameters bounds (do not mind a warning message)
eCopula(rge,ge,technique="ML",method="BFGS")  
# least-square fit to empirical copula
eCopula(rge,ge,technique="LS",procedure="nlminb")  
# maximizing discretized likelihood function
eCopula(rge,ge,technique="ML",procedure="grid",glimits=list(2.,6.),pgrid=20)  
# specify nodes of the grid
eCopula(rge,ge,tech="ML",proc="grid",ggridparameters=list(c(2.,6.,length.out=20))) 
# without naming, it won't create sequence
eCopula(rge,ge,technique="ML",procedure="grid",ggridparameters=list(c(2.,6.,20)))
# inversion of Kendall's tau
eCopula(rge,ge,technique="icorr",corrtype="kendall")

## --- GoF test, set higher N to increase precision of p-value
gCopula(rge,ge,etechnique="ML",N=10)
# parallel computing takes lesser time, but the progress is not displayed
# not available on Windows OS
if(.Platform$OS.type!="windows") {
  gCopula(rge,ge,etechnique="ML",N=10,ncores=2)
}

## testing if two data sets has equal copulas
rge1 <- rCopula(80,dim=2,genClayton(),gpars=3)
gCopula(list(rge,rge1),N=10)

## check whether some hypotheticaly-copula function does not violate  
## copula properties (over data and parameters grid)
isCopula(genGumbel(),dagrid=10,pgrid=10,tolerance=1e-15)

## all the above functions are ready for archimax or generic copulas 
## as well as for higher dimensions
pCopula(c(0.3,0.5,1.0),genClayton(),depGumbel(),gpars=0.01,dpars=4.)  #0.2907613
pCopula(c(0.3,0.5,1.0),copula=copGumbel(),pars=4.)  #0.2906142




cleanEx()
nameEx("xPareto")
### * xPareto

flush(stderr()); flush(stdout())

### Name: xPareto
### Title: 4-parametric univariate Pareto distribution
### Aliases: pPareto qPareto
### Keywords: Pareto distribution CDF quantile

### ** Examples

## probability P(X<q)=p
pPareto(t = 2.5, pars = c(10.,5.,3.,1))  # 0.8823436
qPareto(t = .Last.value, pars = c(10.,5.,3.,1))  # 2.5



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
