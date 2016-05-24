###############################################################################
## Example: Gamma Family
###############################################################################
require(ROptEst)
options("newDevice"=TRUE)


## generates Gamma Family with 
## scale = 2 and shape = 0.1
G <- GammaFamily(scale = 2, shape = 0.1)
G       # show G
plot(G) # plot of Gammad(scale = 2, shape = 0.1) and L_2 derivative
distrExOptions(ErelativeTolerance = 1e-8) # increase precision for E
checkL2deriv(G)

## 30.06.09: new method for "E" with 
## signature(object = "Gammad", fun = "function", cond = "missing")
## in package distrEx introduced which slightly reduces the problem
## documented below.

## more precisely:
## numerical integration gives
E(Gammad(scale = 2, shape = 0.1), function(x) (log(x/2)-digamma(0.1))^2)

## whereas
trigamma(0.1)

## Problem is more or less caused by integration of log(x/2)^2
## occurs for shape parameter small (<< 1)
E(Gammad(scale = 2, shape = 0.1), function(x) log(x/2)^2)
distrExOptions(ErelativeTolerance = 1e-10)
E(Gammad(scale = 2, shape = 0.1), function(x) log(x/2)^2)
distrExOptions(ErelativeTolerance = 1e-7)
E(Gammad(scale = 2, shape = 0.1), function(x) log(x/2)^2)

## result should be
res1 <- 2*digamma(0.1)*E(Gammad(scale = 2, shape = 0.1), function(x) log(x/2))
res2 <- digamma(0.1)^2
trigamma(0.1) + res1 - res2

fun <- function(x) log(x/2)^2*dgamma(x, scale = 2, shape = 0.1)
integrate(fun, lower = 0, upper = Inf, rel.tol = 1e-6)
integrate(fun, lower = 0, upper = Inf, rel.tol = 1e-3)
integrate(fun, lower = 1e-9, upper = Inf, rel.tol = 1e-6)
integrate(fun, lower = 1e-13, upper = Inf, rel.tol = 1e-6)
## problem at zero!
fun(0)
curve(fun, from = 0, to = 1, n = 501)

fun <- function(x) (log(x/2)-digamma(0.1))^2*dgamma(x, scale = 2, shape = 0.1)
curve(fun, from = 0, to = 1, n = 501)

GLIntegrate(fun, lower = 1e-9, upper = qgamma(1-1e-9, scale = 2, shape = 0.1), order = 500)


## generates Gamma Family with 
## scale = 2 and shape = 1
G <- GammaFamily(scale = 2, shape = 1)
G       # show G
plot(G) # plot of Gammad(scale = 2, shape = 1) and L_2 derivative
distrExOptions(ErelativeTolerance = 1e-8) # increase precision for E
checkL2deriv(G)


## classical optimal IC
IC0 <- optIC(model = G, risk = asCov())
IC0       # show IC
checkIC(IC0) # seems to be a numerical problem!?
Risks(IC0)
plot(IC0) # plot IC

## L_2 family + infinitesimal neighborhood
RobG1 <- InfRobModel(center = G, neighbor = ContNeighborhood(radius = 0.5))
RobG1     # show RobB1

## OBRE solution ARE 0.95 in ideal model
system.time(ICA <- optIC(model = RobG1, risk = asAnscombe(),
                         upper=NULL,lower=NULL, verbose=TRUE))
checkIC(ICA)
Risks(ICA)
plot(ICA)
infoPlot(ICA)


## MSE solution
system.time(IC1 <- optIC(model=RobG1, risk=asMSE()))
IC1
checkIC(IC1)
Risks(IC1)
plot(IC1)
#devNew()
infoPlot(IC1)

## lower case solutions
system.time(IC2 <- optIC(model=RobG1, risk=asBias(), tol = 1e-10))
IC2
checkIC(IC2)
Risks(IC2)
plot(IC2)
devNew()
infoPlot(IC2)

## Hampel solution
system.time(IC3 <- optIC(model=RobG1, risk=asHampel(bound=clip(IC1))))
IC3
checkIC(IC3)
Risks(IC3)
plot(IC3)
x11()
infoPlot(IC3)

## radius minimax IC
## takes quite some time - about 30 min.
system.time(IC4 <- radiusMinimaxIC(L2Fam=G, neighbor=ContNeighborhood(), 
            risk=asMSE(), loRad=0, upRad=Inf))

## least favorable radius
## takes really long time - several hours!
#system.time(r.rho1 <- leastFavorableRadius(L2Fam=G, neighbor=ContNeighborhood(),
#                    risk=asMSE(), rho=0.5))

## one-step estimation
## 1. generate a contaminated sample
ind <- rbinom(100, size=1, prob=0.05) 
x <- (1-ind)*rgamma(100, scale = 1, shape = 2) + ind*rgamma(100, scale = 3, shape = 5)

## 2. Kolmogorov(-Smirnov) minimum distance estimator
(est01 <- MDEstimator(x=x, GammaFamily()))

## 3. Cramer von Mises minimum distance estimator
(est02 <- MDEstimator(x = x, GammaFamily(), distance = CvMDist))

## non-robust ML estimator
MLEstimator(x=x, GammaFamily())

## 4. one-step estimation: radius known
RobG3 <- InfRobModel(center=GammaFamily(scale = estimate(est02)["scale"], 
                     shape = estimate(est02)["shape"]), 
                neighbor=ContNeighborhood(radius=0.5))
IC9 <- optIC(model=RobG3, risk=asMSE())
(est1 <- oneStepEstimator(x, IC=IC9, start=est02))

## 5. k-step estimation: radius known
(est2 <- kStepEstimator(x, IC=IC9, start=est02, steps = 3))

## It's simpler to use function roptest!
(est3 <- roptest(x, eps = 0.5/sqrt(length(x)), L2Fam = GammaFamily()))
(est4 <- roptest(x, eps = 0.5/sqrt(length(x)), L2Fam = GammaFamily(), steps = 3))

## comparison
estimate(est1)
estimate(est3)
estimate(est2)
estimate(est4)

## confidence intervals
confint(est1, symmetricBias())
confint(est3, symmetricBias())
confint(est2, symmetricBias())
confint(est4, symmetricBias())
