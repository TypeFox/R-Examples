### R code from vignette source 'distr.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: SweaveListingsPreparations
###################################################
require(SweaveListingUtils)
SweaveListingPreparations()
setToBeDefinedPkgs(pkgs = c("distr","distrEx","distrTEst","distrSim",
                            "distrDoc","distrTeach","distrMod","RandVar"),
                   keywordstyles = "\\bf\\color{distrCol}")


###################################################
### code chunk number 2: distr.Rnw:229-236
###################################################
## preparation: set option withSweave to true
require(distrTEst)
require(distrEx)
require(distrTeach)
require(distrMod)
distroptions(withSweave = TRUE)
options("newDevice" = TRUE)


###################################################
### code chunk number 3: exam1
###################################################
require(distr)
N <- Norm(mean = 2, sd = 1.3)
P <- Pois(lambda = 1.2)
Z <- 2*N + 3 + P
Z
plot(Z)
p(Z)(0.4)
q(Z)(0.3)
Zs <- r(Z)(50)
Zs


###################################################
### code chunk number 4: DiscrDist
###################################################
D <- DiscreteDistribution(supp = c(1,5,7,21), prob = c(0.1,0.1,0.6,0.2))
D
plot(D)


###################################################
### code chunk number 5: AbscDist
###################################################
AC <- AbscontDistribution(d = function(x) exp(-abs(x)^3), withStand = TRUE)
AC
plot(AC)


###################################################
### code chunk number 6: examLis
###################################################
library(distr)
M1 <- UnivarMixingDistribution(Norm(), Pois(lambda=1), Norm(), 
      withSimplify = FALSE)
M2 <- UnivarMixingDistribution(M1, Norm(), M1, Norm(), withSimplify = FALSE)
M2


###################################################
### code chunk number 7: warningArith
###################################################
  A1 <- Norm(); A2 <- Unif()
  A1 + A2


###################################################
### code chunk number 8: examdcP
###################################################
decomposePM(Norm())
     decomposePM(Binom(2,0.3)-Binom(5,.4))
     decomposePM(UnivarLebDecDistribution(Norm(),Binom(2,0.3)-Binom(5,.4), 
                 acWeight = 0.3))


###################################################
### code chunk number 9: examflat
###################################################
D1 <- Norm()
D2 <- Pois(1)
D3 <- Binom(1,.4)
D4 <- UnivarMixingDistribution(D1,D2,D3, mixCoeff = c(0.4,0.5,0.1), 
      withSimplify = FALSE)
D <- UnivarMixingDistribution(D1,D4,D1,D2, mixCoeff = c(0.4,0.3,0.1,0.2), 
      withSimplify = FALSE)
D
D0<-flat.mix(D)
D0


###################################################
### code chunk number 10: arithmetic
###################################################
  A1 <- Norm(); A2 <- Unif()
  d(sin(A1 + A2))(0.1)
  d(sin(A1 + A2))(0.1)
  sin(A1 + A2)


###################################################
### code chunk number 11: arith2v1
###################################################
  A1 <- Norm(); A2 <- Unif()
  A1A2 <- A1*A2
  plot(A1A2)


###################################################
### code chunk number 12: arith2v2
###################################################
  A12 <- 1/(A2 + .3)
  plot(A12) 


###################################################
### code chunk number 13: arith2v3
###################################################
  B <- Binom(5,.2)+1
  A1B <- A1^B
  plot(A1B, xlim=c(-3,3))


###################################################
### code chunk number 14: arith2V4
###################################################
  plot(1.2^A1)


###################################################
### code chunk number 15: arith2v5
###################################################
  plot(B^A1)


###################################################
### code chunk number 16: Hub
###################################################
H <- Huberize(Norm(),lower=-1,upper=2)
plot(H)


###################################################
### code chunk number 17: Trun
###################################################
T <- Truncate(Norm(),lower=-1,upper=2)
plot(T)


###################################################
### code chunk number 18: Min1
###################################################
M1 <- Maximum(Unif(0,1), Minimum(Unif(0,1), Unif(0,1)))
plot(M1)


###################################################
### code chunk number 19: Min2
###################################################
M2 <- Minimum(Exp(4),4)
plot(M2)


###################################################
### code chunk number 20: Min3
###################################################
M3 <- Minimum(Norm(2,2), Pois(3))
plot(M3)


###################################################
### code chunk number 21: TruncExtr
###################################################
N <- Norm()
TN <- Truncate(N, 20,22)
r(TN)(20)  ## simulation accurate although :
p(N)(20, lower.tail = FALSE) ## prob that N>=20 


###################################################
### code chunk number 22: qrex
###################################################
B <- Binom(5,0.5)
p(B)(3)
p.l(B)(3)
q(B)(.5)
q.r(B)(0.5)


###################################################
### code chunk number 23: probHN
###################################################
B0 <- as(Binom(5,0.5),"DiscreteDistribution")
   ## coercion necessary:
   ## otherwise slot "prob" of B0 will be returned
prob(B0)
HN <- Huberize(N, -2,1)
prob(HN)


###################################################
### code chunk number 24: makeAC
###################################################
par(mfrow=c(2,3))
plot(makeAbscontDistribution(Nbinom(5,.5)),mfColRow=FALSE)
plot(makeAbscontDistribution(HN),mfColRow=FALSE)
par(mfrow=c(1,1))


###################################################
### code chunk number 25: getLowUp
###################################################
getLow(Nbinom(5,0.5))
getUp(Nbinom(5,0.5))
getLow(Norm(5,0.5))
getUp(Norm(5,0.5))


###################################################
### code chunk number 26: cauchy1
###################################################
  plot(Cauchy())


###################################################
### code chunk number 27: cauchy2
###################################################
  plot(Cauchy(),xlim=c(-4,4))


###################################################
### code chunk number 28: plotex1
###################################################
plot(Binom(size = 4, prob = 0.3))


###################################################
### code chunk number 29: plotex2
###################################################
plot(Binom(size = 4, prob = 0.3), do.points = FALSE, verticals = FALSE)


###################################################
### code chunk number 30: plotex3
###################################################
plot(Binom(size = 4, prob = 0.3), main = TRUE, inner = FALSE, cex.main = 1.6,
     tmar = 6)


###################################################
### code chunk number 31: plotex4
###################################################
plot(Binom(size = 4, prob = 0.3), cex.points = 1.2, pch = 20, lwd = 2)


###################################################
### code chunk number 32: plotex5
###################################################
B <- Binom(size = 4, prob = 0.3)
plot(B, col="red", col.points = "green", main = TRUE, col.main="blue",
     col.sub = "orange", sub = TRUE, cex.sub = 0.6, col.inner = "brown")


###################################################
### code chunk number 33: plotex6
###################################################
plot(Nbinom(size = 4,prob = 0.3), cex.points = 1.2, pch.u = 20, pch.a = 10)


###################################################
### code chunk number 34: plotex7
###################################################
plot(Chisq(), log = "xy", ngrid = 100)


###################################################
### code chunk number 35: plotex8
###################################################
plot(Norm(), lwd=3, col = "red", ngrid = 200, lty = 3, las = 2)


###################################################
### code chunk number 36: plotex9
###################################################
plot(Norm(), panel.first = grid(), main = "my Distribution: %A",
     inner = list(expression(paste(lambda, "-density of %C(%P)")), "CDF",
                  "Pseudo-inverse with param's %N"),
     sub = "this plot was correctly generated on %D",
     cex.inner = 0.9, cex.sub = 0.8)


###################################################
### code chunk number 37: plotex10
###################################################
Ch <- Chisq(); setgaps(Ch, exactq = 3)
plot(Ch, cex = 1.2, pch.u = 20, pch.a = 10, col.points = "green", 
     col.vert = "red")


###################################################
### code chunk number 38: plotex11
###################################################
layout(matrix(c(1,3,2,3), nrow=2))
plot(N, mfColRow = FALSE)


###################################################
### code chunk number 39: plotex12
###################################################
layout(matrix(c(rep(1,6),2,2,3,3,4,4,5,5,5,6,6,6), 
                   nrow=3, byrow=TRUE))
plot(HN, mfColRow = FALSE,
        to.draw.arg=c("p","d.c","p.c","q.c", "p.d","q.d"))


###################################################
### code chunk number 40: simulate
###################################################
X <- Simulation()
seed(X) <- setRNG()
simulate(X)
Data(X)[1:10]


###################################################
### code chunk number 41: expectation
###################################################
D4 <- LMCondDistribution(theta = 1)
D4  # corresponds to Norm(cond, 1)
N <- Norm(mean = 2)

E(D4, cond = 1)
E(D4, cond = 1, useApply = FALSE)
E(as(D4, "UnivariateCondDistribution"), cond = 1)
E(D4, function(x){x^2}, cond = 2)
E(D4, function(x){x^2}, cond = 2, useApply = FALSE)
E(N, function(x){x^2})
E(as(N, "UnivariateDistribution"), function(x){x^2}, 
     useApply = FALSE) # crude Monte-Carlo
E(D4, function(x, cond){cond*x^2}, cond = 2,
  withCond = TRUE)
E(D4, function(x, cond){cond*x^2}, cond = 2,
  withCond = TRUE, useApply = FALSE)
E(N, function(x){2*x^2})
E(as(N, "UnivariateDistribution"), function(x){2*x^2},
  useApply = FALSE) # crude Monte-Carlo
Y <- 5 * Binom(4, .25) - 3
Y
E(Y)  


###################################################
### code chunk number 42: expectationlow
###################################################
E(Cauchy(), low=3, upp=5)
var(Cauchy(), low=3, upp=5)


###################################################
### code chunk number 43: expectation2
###################################################
E(N, function(x)x^2) 
E(N, function(x)x^2,  lowerTruncQuantile = 1e-5)
var(Cauchy(), low =3, upperTruncQuantile = 1e-5,  IQR.fac = 10)
var(Cauchy(), low =3, upperTruncQuantile = 1e-10, IQR.fac = 20)


###################################################
### code chunk number 44: var (eval = FALSE)
###################################################
##   var <- function(x , ...)
##        {dots <- list(...)
##         if(hasArg(y)) y <- dots$"y"
##         na.rm <- ifelse(hasArg(na.rm), dots$"na.rm", FALSE)
##         if(!hasArg(use))
##              use <- ifelse (na.rm, "complete.obs","all.obs")
##         else use <- dots$"use"
##         if(hasArg(y))
##            stats::var(x = x, y = y, na.rm = na.rm, use)
##         else
##            stats::var(x = x, y = NULL, na.rm = na.rm, use)
##         }


###################################################
### code chunk number 45: MCEstimator
###################################################
    library(distrMod)
    x <- rgamma(50, scale = 0.5, shape = 3)
    G <- GammaFamily(scale = 1, shape = 2)
    negLoglikelihood <- function(x, Distribution){
        res <- -sum(log(Distribution@d(x)))
        names(res) <- "Negative Log-Likelihood"
        return(res)
    }
    MCEstimator(x = x, ParamFamily = G, criterion = negLoglikelihood)


###################################################
### code chunk number 46: censPoisFamilyDef (eval = FALSE)
###################################################
##     ## search interval for reasonable parameters
##     startPar <- function(x,...) c(.Machine$double.eps,max(x))
## 
##     ## what to do in case of leaving the parameter domain
##     makeOKPar <- function(param) {if(param<=0) return(.Machine$double.eps)
##                                   return(param)}


###################################################
### code chunk number 47: PoisFamilyDef (eval = FALSE)
###################################################
## setClass("PoisFamily", contains = "L2ParamFamily")


###################################################
### code chunk number 48: NormLocationFamily (eval = FALSE)
###################################################
## setClass("NormLocationFamily", contains = "L2LocationFamily")


###################################################
### code chunk number 49: L2ScaleFamily (eval = FALSE)
###################################################
##  setMethod("validParameter", signature(object = "L2ScaleFamily"),
##           function(object, param, tol=.Machine$double.eps){
##              if(is(param,"ParamFamParameter"))
##                 param <- main(param)
##              if(!all(is.finite(param))) return(FALSE)
##              if(length(param)!=1) return(FALSE)
##              return(param > tol)})


###################################################
### code chunk number 50: GammaFamilyModify (eval = FALSE)
###################################################
## setMethod("modifyModel", signature(model = "GammaFamily",
##            param = "ParamFamParameter"),
##           function(model, param, ...){
##              M <- modifyModel(as(model, "L2ParamFamily"), param = param,
##                               .withCall = FALSE)
##              M@L2derivSymm <- FunSymmList(OddSymmetric(SymmCenter = 
##                                                        prod(main(param))),
##                                           NonSymmetric())
##              class(M) <- class(model)
##              return(M)
##           })


###################################################
### code chunk number 51: MLEstimator
###################################################
    MLEstimator(x = x, ParamFamily = G)
    MDEstimator(x = x, ParamFamily = G, distance = KolmogorovDist)


###################################################
### code chunk number 52: NormScaleFam (eval = FALSE)
###################################################
## setMethod("mleCalc", signature(x = "numeric", PFam = "NormScaleFamily"),
##            function(x, PFam, ...){
##            n <- length(x)
##            theta <- sqrt((n-1)/n)*sd(x); mn <- mean(distribution(PFam))
##            ll <- -sum(dnorm(x, mean=mn, sd = theta, log=TRUE))
##            names(ll) <- "neg.Loglikelihood"
##            crit.fct <- function(sd)
##                          -sum(dnorm(x, mean=mn, sd = sd, log=TRUE))  
##            param <- ParamFamParameter(name = "scale parameter", 
##                                main = c("sd"=theta))
##            if(!hasArg(Infos)) Infos <- NULL
##            return(meRes(x, theta, ll, param, crit.fct, Infos = Infos))
## })


###################################################
### code chunk number 53: CIex
###################################################
require(distrMod)
## some transformation
mtrafo <- function(x){
     nms0 <- c("scale","shape")
     nms <- c("shape","rate")
     fval0 <- c(x[2], 1/x[1])
     names(fval0) <- nms
     mat0 <- matrix( c(0, -1/x[1]^2, 1, 0), nrow = 2, ncol = 2,
                     dimnames = list(nms,nms0))                          
     list(fval = fval0, mat = mat0)}

set.seed(124)
x <- rgamma(50, scale = 0.5, shape = 3)

## parametric family of probability measures
G <- GammaFamily(scale = 1, shape = 2, trafo = mtrafo)
## MLE
res <- MLEstimator(x = x, ParamFamily = G)
print(res, digits = 4, show.details="maximal")
print(res, digits = 4, show.details="medium")
print(res, digits = 4, show.details="minimal")
ci <- confint(res)
print(ci, digits = 4, show.details="maximal")
print(ci, digits = 4, show.details="medium")
print(ci, digits = 4, show.details="minimal")

## some profiling
par(mfrow=c(2,1))
plot(profile(res))


###################################################
### code chunk number 54: NormApprox
###################################################
require(distr)

N <- Norm(0,1)
U <- Unif(0,1)
U2 <- U + U 
U4 <- U2 + U2
U8 <- U4 + U4
U12 <- U4 + U8
NormApprox <- U12 - 6

x <- seq(-4,4,0.001)

opar <- par(no.readonly = TRUE)
par(mfrow = c(2,1))

plot(x, d(NormApprox)(x),
     type = "l",
     xlab = "",
     ylab = "Density",
     main = "Exact and approximated density")
lines(x, d(N)(x),
      col = "red")
legend("topleft",
       legend = c("NormApprox", "Norm(0,1)"),
       fill = c("black", "red"))

plot(x, d(NormApprox)(x) - d(N)(x),
     type = "l",
     xlab = "",
     ylab = "\"black\" - \"red\"",
     col = "darkgreen",
     main = "Error")
lines(c(-4,4), c(0,0))

par(opar)


###################################################
### code chunk number 55: ConvolutionNormalDistr
###################################################
require(distr)

## initialize two normal distributions
A <- Norm(mean=1, sd=2)
B <- Norm(mean=4, sd=3) 

## convolution via addition of moments
AB <- A+B

## casting of A,B as absolutely continuous distributions
## that is, ``forget'' that A,B are normal distributions
A1 <- as(A, "AbscontDistribution")
B1 <- as(B, "AbscontDistribution")

## for higher precision we change the global variable
## "TruncQuantile" from 1e-5 to 1e-8
oldeps <- getdistrOption("TruncQuantile")
eps <- 1e-8
distroptions("TruncQuantile" = eps)
## support of A1+B1 for FFT convolution is
## [q(A1)(TruncQuantile), 
##  q(B1)(TruncQuantile, lower.tail = FALSE)]

## convolution via FFT
AB1 <- A1+B1

#############################
## plots of the results
#############################
par(mfrow=c(1,3))
low <- q(AB)(1e-15)
upp <- q(AB)(1e-15, lower.tail = FALSE)
x <- seq(from = low, to = upp, length = 10000)

## densities
plot(x, d(AB)(x), type = "l", lwd = 5)
lines(x , d(AB1)(x), col = "orange", lwd = 1)
title("Densities")
legend("topleft", legend=c("exact", "FFT"), 
        fill=c("black", "orange"))

## cdfs
plot(x, p(AB)(x), type = "l", lwd = 5)
lines(x , p(AB1)(x), col = "orange", lwd = 1)
title("CDFs")
legend("topleft", legend=c("exact", "FFT"), 
        fill=c("black", "orange"))

## quantile functions
x <- seq(from = eps, to = 1-eps, length = 1000)
plot(x, q(AB)(x), type = "l", lwd = 5)
lines(x , q(AB1)(x), col = "orange", lwd = 1) 
title("Quantile functions")
legend("topleft", legend=c("exact", "FFT"), 
        fill=c("black", "orange"))

## Since the plots of the results show no 
## recognizable differencies, we also compute 
## the total variation distance of the densities 
## and the Kolmogorov distance of the cdfs

## total variation distance of densities
total.var <- function(z, N1, N2){
    0.5*abs(d(N1)(z) - d(N2)(z))
}
dv <- integrate(total.var, lower=-Inf, upper=Inf, rel.tol=1e-8, N1=AB, N2=AB1)
cat("Total variation distance of densities:\t")
print(dv) # 4.25e-07

### meanwhile realized in package "distrEx" 
### as TotalVarDist(N1,N2)

## Kolmogorov distance of cdfs 
## the distance is evaluated on a random grid
z <- r(Unif(Min=low, Max=upp))(1e5)
dk <- max(abs(p(AB)(z)-p(AB1)(z)))
cat("Kolmogorov distance of cdfs:\t", dk, "\n") 
# 2.03e-07

### meanwhile realized in package "distrEx" 
### as KolmogorovDist(N1,N2)

## old distroptions
distroptions("TruncQuantile" = oldeps)



###################################################
### code chunk number 56: ComparisonFFTandRtoDPQ
###################################################
require(distr)

################################
## Comparison 1 - FFT and RtoDPQ 
################################

N1 <- Norm(0,3)
N2 <- Norm(0,4)
rnew1 <- function(n) r(N1)(n) + r(N2)(n) 

X <- N1 + N2 
     # exact formula -> N(0,5)
Y <- N1 + as(N2, "AbscontDistribution") 
     # appoximated with FFT
Z <- new("AbscontDistribution", r = rnew1) 
     # appoximated with RtoDPQ

# density-plot

x <- seq(-15,15,0.01)
plot(x, d(X)(x),
     type = "l",
     lwd = 3,
     xlab = "",
     ylab = "density",
     main = "Comparison 1",     
     col = "black")
lines(x, d(Y)(x),
     col = "yellow")
lines(x, d(Z)(x),
     col = "red")
legend("topleft",
  legend = c("Exact", "FFT-Approximation", 
             "RtoDPQ-Approximation"),
       fill = c("black", "yellow", "red"))
      
############################################
## Comparison 2 - "Exact" Formula and RtoDPQ
############################################

B <- Binom(size = 6, prob = 0.5) * 10
N <- Norm()
rnew2 <- function(n) r(B)(n) + r(N)(n)

Y <- B + N 
     # "exact" formula
Z <- new("AbscontDistribution", r = rnew2) 
     # appoximated with RtoDPQ

# density-plot

x  <- seq(-5,65,0.01)
plot(x, d(Y)(x),
     type = "l",
     xlab = "",
     ylab = "density",
     main = "Comparison 2",
     col = "black")
lines(x, d(Z)(x),
     col = "red")
legend("topleft",
       legend = c("Exact", "RtoDQP-Approximation"),
       fill = c("black", "red"))


###################################################
### code chunk number 57: StationaryRegressorDistr
###################################################
require(distr)

## Approximation of the stationary regressor 
## distribution of an AR(1) process 
##       X_t = phi X_{t-1} + V_t 
## where V_t i.i.d N(0,1) and phi\in(0,1)
## We obtain 
##    X_t = \sum_{j=1}^\infty phi^j V_{t-j}
## i.e., X_t \sim N(0,1/(1-phi^2))
phi <- 0.5

## casting of V as absolutely continuous distributions
## that is, ``forget'' that V is a normal distribution
V <- as(Norm(), "AbscontDistribution")

## for higher precision we change the global variable
## "TruncQuantile" from 1e-5 to 1e-8
oldeps <- getdistrOption("TruncQuantile")
eps <- 1e-8
distroptions("TruncQuantile" = eps)

## Computation of the approximation 
##      H=\sum_{j=1}^n phi^j V_{t-j}
## of the stationary regressor distribution 
## (via convolution using FFT)
H <- V
n <- 15 
## may take some time
### switch off warnings [would be issued due to 
###  very unequal variances...]
old.warn <- getOption("warn")
options("warn" = -1)
for(i in 1:n){Vi <- phi^i*V; H <- H + Vi } 
options("warn" = old.warn)

## the stationary regressor distribution (exact)
X <- Norm(sd=sqrt(1/(1-phi^2)))

#############################
## plots of the results
#############################
par(mfrow=c(1,3))
low <- q(X)(1e-15)
upp <- q(X)(1e-15, lower.tail = FALSE)
x <- seq(from = low, to = upp, length = 10000)

## densities
plot(x, d(X)(x),type = "l", lwd = 5)
lines(x , d(H)(x), col = "orange", lwd = 1)
title("Densities")
legend("topleft", legend=c("exact", "FFT"), 
        fill=c("black", "orange"))

## cdfs
plot(x, p(X)(x),type = "l", lwd = 5)
lines(x , p(H)(x), col = "orange", lwd = 1)
title("CDFs")
legend("topleft", legend=c("exact", "FFT"), 
        fill=c("black", "orange"))

## quantile functions
x <- seq(from = eps, to = 1-eps, length = 1000)
plot(x, q(X)(x),type = "l", lwd = 5)
lines(x , q(H)(x), col = "orange", lwd = 1)
title("Quantile functions")
legend( "topleft", 
        legend=c("exact", "FFT"), 
        fill=c("black", "orange"))

## Since the plots of the results show no 
## recognizable differencies, we also compute 
## the total variation distance of the densities 
## and the Kolmogorov distance of the cdfs

## total variation distance of densities
total.var <- function(z, N1, N2){
    0.5*abs(d(N1)(z) - d(N2)(z))
}
dv <- integrate(f = total.var, lower = -Inf, 
                upper = Inf, rel.tol = 1e-7, 
                N1=X, N2=H)
cat("Total variation distance of densities:\t")
print(dv) # ~ 5.0e-06

### meanwhile realized in package "distrEx" 
### as TotalVarDist(N1,N2)


## Kolmogorov distance of cdfs 
## the distance is evaluated on a random grid
z <- r(Unif(Min=low, Max=upp))(1e5)
dk <- max(abs(p(X)(z)-p(H)(z)))
cat("Kolmogorov distance of cdfs:\t", dk, "\n") 
# ~2.5e-06

### meanwhile realized in package "distrEx" 
### as KolmogorovDist(N1,N2)


## old distroptions
distroptions("TruncQuantile" = oldeps)



###################################################
### code chunk number 58: destructive
###################################################
##########################################################
## Demo: Instructive destructive example
##########################################################
require(distr)

## package "distr" encourages 
## consistency but does not 
## enforce it---so in general  
## d o   n o t   m o d i f y
## slots d,p,q,r!

N <- Norm()
B <- Binom()
N@d <- B@d
plot(N, lwd = 3) 


###################################################
### code chunk number 59: SimulateandEstimate
###################################################
require(distrTEst)
    ### also loads distrSim
sim <- new("Simulation",
           seed = setRNG(),
           distribution = Norm(mean = 0, sd = 1),
           filename="sim_01",
           runs = 1000,
           samplesize = 30)

contsim <- new("Contsimulation",
               seed = setRNG(),
               distribution.id = Norm(mean = 0, sd = 1),
               distribution.c = Norm(mean = 0, sd = 9),
               rate = 0.1,
               filename="contsim_01",
               runs = 1000,
               samplesize = 30)

simulate(sim)
simulate(contsim)

sim
summary(contsim)
plot(contsim)


###################################################
### code chunk number 60: elist
###################################################
require(distrTEst)
psim <- function(theta,y,m0){
  mean(pmin(pmax(-m0, y - theta), m0))
  }
mestimator <- function(x, m = 0.7) {
  uniroot(f = psim,
          lower = -20,
          upper = 20,
          tol = 1e-10,
          y = x,
          m0 = m,
          maxiter = 20)$root
  }

result.id.mean <- evaluate(sim, mean)
result.id.mest <- evaluate(sim, mestimator)
result.id.median <- evaluate(sim, median)


result.cont.mean <- evaluate(contsim, mean)
result.cont.mest <- evaluate(contsim, mestimator)
result.cont.median <- evaluate(contsim, median)

elist <- EvaluationList(result.cont.mean,
                        result.cont.mest,
                        result.cont.median) 

elist
summary(elist)
plot(elist, cex = 0.7, las = 2)


###################################################
### code chunk number 61: mest
###################################################
result.cont.mest


###################################################
### code chunk number 62: elist.summary
###################################################
summary(elist)


###################################################
### code chunk number 63: Expectation
###################################################
require("distrEx")
# Example
id <- function(x) x
sq <- function(x) x^2

# Expectation and Variance of Binom(6,0.5)
B <- Binom(6, 0.5)
E(B, id)
E(B, sq) - E(B, id)^2

# Expectation and Variance of Norm(1,1)
N <- Norm(1, 1)
E(N, id)
E(N, sq) - E(N, id)^2


###################################################
### code chunk number 64: nFoldConvolution
###################################################
##########################################################
## Demo: n-fold convolution of absolutely continuous 
##       probability distributions
##########################################################
require(distr)

if(!isGeneric("convpow")) 
    setGeneric("convpow", 
    function(D1,...) standardGeneric("convpow"))

##########################################################
## Function for n-fold convolution
## -- absolute continuous distribution -- 
##########################################################

##implentation of Algorithm 3.4. of 
# Kohl, M., Ruckdeschel, P., Stabla, T. (2005): 
#   General purpose convolution algorithm for distributions 
#   in S4-Classes by means of FFT.
# Technical report, Feb. 2005. Also available in
# http://www.uni-bayreuth.de/departments/math/org/mathe7/
#       /RUCKDESCHEL/pubs/comp.pdf


setMethod("convpow",
          signature(D1 = "AbscontDistribution"),
          function(D1, N){
            if((N < 1)||(!identical(floor(N), N)))
              stop("N has to be a natural greater than 0")
            
            m <- getdistrOption("DefaultNrFFTGridPointsExponent")

    ##STEP 1

            lower <- ifelse((q(D1)(0) > - Inf), q(D1)(0), 
                     q(D1)(getdistrOption("TruncQuantile"))) 
            upper <- ifelse((q(D1)(1) < Inf), q(D1)(1), 
                     q(D1)(getdistrOption("TruncQuantile"), lower.tail = FALSE))

    ##STEP 2

            M <- 2^m
            h <- (upper-lower)/M
            if(h > 0.01)
              warning(paste("Grid for approxfun too wide, ", 
              "increase DefaultNrFFTGridPointsExponent", sep=""))
            x <- seq(from = lower, to = upper, by = h)
            p1 <- p(D1)(x)

    ##STEP 3

            p1 <- p1[2:(M + 1)] - p1[1:M]

    ##STEP 4
    
            ## computation of DFT
            pn <- c(p1, numeric((N-1)*M))
            fftpn <- fft(pn)

    ##STEP 5

            ## convolution theorem for DFTs
            pn <- Re(fft(fftpn^N, inverse = TRUE)) / (N*M)
            pn <- (abs(pn) >= .Machine$double.eps)*pn
            i.max <- N*M-(N-2)
            pn <- c(0,pn[1:i.max])
            dn <- pn / h
            pn <- cumsum(pn)

    ##STEP 6(density)

            ## density 
            x <- c(N*lower,seq(from = N*lower+N/2*h, 
                   to = N*upper-N/2*h, by=h),N*upper)
            dnfun1 <- approxfun(x = x, y = dn, yleft = 0, yright = 0)

    ##STEP 7(density)
 
            standardizer <- sum(dn[2:i.max]) + (dn[1]+dn[i.max+1]) / 2
            dnfun2 <- function(x) dnfun1(x) / standardizer

    ##STEP 6(cdf)
    
            ## cdf with continuity correction h/2
            pnfun1 <- approxfun(x = x+0.5*h, y = pn, 
                        yleft = 0, yright = pn[i.max+1])

    ##STEP 7(cdf)
   
            pnfun2 <- function(x) pnfun1(x) / pn[i.max+1]


            ## quantile with continuity correction h/2
            yleft <- ifelse(((q(D1)(0) == -Inf)|
                             (q(D1)(0) == -Inf)), 
                             -Inf, N*lower)
            yright <- ifelse(((q(D1)(1) == Inf)|
                              (q(D1)(1) == Inf)), 
                              Inf, N*upper)    
            w0 <- options("warn")
            options(warn = -1)
            qnfun1 <- approxfun(x = pnfun2(x+0.5*h), 
                        y = x+0.5*h, yleft = yleft, yright = yright)
            qnfun2 <- function(x){ 
            ind1 <- (x == 0)*(1:length(x))
            ind2 <- (x == 1)*(1:length(x))
            y <- qnfun1(x)
            y <- replace(y, ind1[ind1 != 0], yleft)
            y <- replace(y, ind2[ind2 != 0], yright)
            return(y)
            }
            options(w0)

            rnew = function(N) apply(matrix(r(e1)(n*N), 
                                     ncol=N), 1, sum)

            return(new("AbscontDistribution", r = rnew, 
                       d = dnfun1, p = pnfun2, q = qnfun2))
})


## initialize a normal distribution
A <- Norm(mean=0, sd=1)

## convolution power
N <- 10 

## convolution via FFT
AN <- convpow(as(A,"AbscontDistribution"), N)
##  ... for the normal distribution , 'convpow' has an "exact"
##      method by version 1.9 so the as(.,.)  is needed to
##      see how the algorithm above works

## convolution exact
AN1 <- Norm(mean=0, sd=sqrt(N))

## plots of the results
eps <- getdistrOption("TruncQuantile")
par(mfrow=c(1,3))
low <- q(AN1)(eps)
upp <- q(AN1)(eps, lower.tail = FALSE)
x <- seq(from = low, to = upp, length = 10000)

## densities
plot(x, d(AN1)(x), type = "l", lwd = 5)
lines(x , d(AN)(x), col = "orange", lwd = 1)
title("Densities")
legend("topleft", legend=c("exact", "FFT"), 
        fill=c("black", "orange"))

## cdfs
plot(x, p(AN1)(x), type = "l", lwd = 5)
lines(x , p(AN)(x), col = "orange", lwd = 1)
title("CDFs")
legend("topleft", legend=c("exact", "FFT"), 
        fill=c("black", "orange"))

## quantile functions
x <- seq(from = eps, to = 1-eps, length = 1000)
plot(x, q(AN1)(x), type = "l", lwd = 5)
lines(x , q(AN)(x), col = "orange", lwd = 1) 
title("Quantile functions")
legend("topleft", 
       legend = c("exact", "FFT"), 
        fill = c("black", "orange"))


