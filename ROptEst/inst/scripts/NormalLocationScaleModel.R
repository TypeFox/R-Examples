###############################################################################
## Example: Normal location and scale
###############################################################################
require(ROptEst)
options("newDevice"=TRUE)

## generates normal location and scale family with mean = -2 and sd = 3
N0 <- NormLocationScaleFamily(mean=-2, sd=3) 
N0        # show G0
plot(N0)  # plot of Norm(mean = -2, sd = 3) and L_2 derivative
checkL2deriv(N0)

## classical optimal IC
N0.IC0 <- optIC(model = N0, risk = asCov())
N0.IC0       # show IC
checkIC(N0.IC0)
Risks(N0.IC0)
plot(N0.IC0) # plot IC

## L_2 family + infinitesimal neighborhood
N0.Rob1 <- InfRobModel(center = N0, neighbor = ContNeighborhood(radius = 0.5))
N0.Rob1     # show N0.Rob1

## OBRE solution (ARE = .95)
system.time(N0.ICA <- optIC(model = N0.Rob1, risk = asAnscombe(), upper=NULL,lower=NULL, verbose=TRUE))
checkIC(N0.ICA)
Risks(N0.ICA)
plot(N0.ICA)
infoPlot(N0.ICA)

system.time(N0.ICA.i <- optIC(model = N0.Rob1, risk = asAnscombe(eff=0.95, normtype=InfoNorm()), upper=NULL,lower=NULL, verbose=TRUE))

## MSE solution
system.time(N0.IC1 <- optIC(model = N0.Rob1, risk = asMSE()))
checkIC(N0.IC1)
Risks(N0.IC1)
plot(N0.IC1)
infoPlot(N0.IC1)

system.time(N0.IC1.i <- optIC(model = N0.Rob1, risk = asMSE(normtype=InfoNorm())))
checkIC(N0.IC1.i)
Risks(N0.IC1.i)
plot(N0.IC1.i)
infoPlot(N0.IC1.i)

system.time(N0.IC1.s <- optIC(model = N0.Rob1, risk = asMSE(normtype=SelfNorm())))
checkIC(N0.IC1.s)
Risks(N0.IC1.s)
plot(N0.IC1.s)
infoPlot(N0.IC1.s)

comparePlot(N0.IC1,N0.IC1.i,N0.IC1.s)


## lower case solutions
(N0.IC2 <- optIC(model = N0.Rob1, risk = asBias(), tol = 1e-10))
checkIC(N0.IC2)
Risks(N0.IC2)
plot(N0.IC2)
infoPlot(N0.IC2)

(N0.IC2.i <- optIC(model = N0.Rob1, risk = asBias(normtype=InfoNorm()), tol = 1e-10))
checkIC(N0.IC2.i)
Risks(N0.IC2.i)
plot(N0.IC2.i)
infoPlot(N0.IC2.i)

## Hampel solution
(N0.IC3 <- optIC(model = N0.Rob1, risk = asHampel(bound = clip(N0.IC1))))
checkIC(N0.IC3)
Risks(N0.IC3)
plot(N0.IC3) 
infoPlot(N0.IC3)

## radius minimax IC
## (may take quite some time!)
system.time(N0.IC4 <- radiusMinimaxIC(L2Fam=N0, neighbor=ContNeighborhood(), 
                risk=asMSE(), loRad=0, upRad=Inf, verbose = TRUE))
checkIC(N0.IC4)
Risks(N0.IC4)
plot(N0.IC4) 
infoPlot(N0.IC4)

(N0.IC4.i <- radiusMinimaxIC(L2Fam=N0, neighbor=ContNeighborhood(), 
                risk=asMSE(normtype=InfoNorm()), loRad=0, upRad=Inf))
checkIC(N0.IC4.i)
Risks(N0.IC4.i)
plot(N0.IC4.i) 
infoPlot(N0.IC4.i)

## takes extremely long time:
(N0.IC4.s <- radiusMinimaxIC(L2Fam=N0, neighbor=ContNeighborhood(), 
                risk=asMSE(normtype=SelfNorm()), loRad=0, upRad=Inf))
checkIC(N0.IC4.s)
Risks(N0.IC4.s)
plot(N0.IC4.s) 
infoPlot(N0.IC4.s)


## least favorable radius
## (may take quite some time!)
N0.r.rho1 <- leastFavorableRadius(L2Fam=N0, neighbor=ContNeighborhood(),
                    risk=asMSE(), rho=0.5)
###############################################################################
# a non-trivial trafo:
###############################################################################

tfct <- function(x){
    nms0 <- c("mean","sd")
    nms  <- "comb"
    fval0 <- x[1]+2*x[2]
    names(fval0) <- nms
    mat0 <- matrix(c(1,2), nrow = 1, dimnames = list(nms,nms0))
    return(list(fval = fval0, mat = mat0))
}

## corresponding ideal and robust models
N1.traf <- NormLocationScaleFamily(mean = 0, sd = 1, trafo= tfct)
N1R.traf <- InfRobModel(center = N1.traf, neighbor = ContNeighborhood(radius = 0.5))
N2R.traf <- InfRobModel(center = N1.traf, neighbor = TotalVarNeighborhood(radius = 0.5))

### classical solution
IC.traf.class <- optIC(model=N1.traf,risk=asCov())
plot(IC.traf.class)
checkIC(IC.traf.class)

### Hampel solution *=c
IC.traf.CV.H <- optIC(model = N1R.traf, risk = asHampel(bound = 8),verbose=TRUE)
plot(IC.traf.CV.H)
checkIC(IC.traf.CV.H)

### MSE solution *=c
IC.traf.CV.MSE <- optIC(model = N1R.traf, risk = asMSE(),verbose=TRUE)
plot(IC.traf.CV.MSE)
checkIC(IC.traf.CV.MSE)

### lower case solution *=c
IC.traf.CV.BIAS <- optIC(model = N1R.traf, risk = asBias(),verbose=TRUE)
plot(IC.traf.CV.BIAS)
checkIC(IC.traf.CV.BIAS)

### Hampel solution *=v
IC.traf.TV.H <- optIC(model = N2R.traf, risk = asHampel(bound = 6),
                      verbose=TRUE, checkBounds=FALSE)
plot(IC.traf.TV.H)
checkIC(IC.traf.TV.H)

### MSE solution *=v
IC.traf.TV.MSE <- optIC(model = N2R.traf, risk = asMSE(),verbose=TRUE)
plot(IC.traf.TV.MSE)
checkIC(IC.traf.TV.MSE)

### lower case solution *=v
IC.traf.TV.BIAS <- optIC(model = N2R.traf, risk = asBias(),verbose=TRUE)
plot(IC.traf.TV.BIAS)
checkIC(IC.traf.TV.BIAS)


###############################################################################
## one-step estimation
###############################################################################
## 1. generate a contaminated sample
ind <- rbinom(100, size=1, prob=0.05) 
x <- rnorm(100, mean=0, sd=(1-ind) + ind*9)

## 2. Kolmogorov(-Smirnov) minimum distance estimator
(est0 <- MDEstimator(x=x, NormLocationScaleFamily()))

## 3. k-step estimation: radius known
N1 <- NormLocationScaleFamily(mean=estimate(est0)["mean"], sd=estimate(est0)["sd"])
N1.Rob <- InfRobModel(center = N1, neighbor = ContNeighborhood(radius = 0.5))
IC1 <- optIC(model = N1.Rob, risk = asMSE())
(est1 <- kStepEstimator(x, IC1, est0, steps = 3))

## 4. k-step estimation: radius unknown
## rough estimate: 1-10% contamination
## => r\in[0.1,1.0]

## takes some time
IC2 <- radiusMinimaxIC(L2Fam=N1, neighbor=ContNeighborhood(),risk=asMSE(), 
                       loRad=0.1, upRad=1.0) 
(est2 <- oneStepEstimator(x, IC2, est0))

########### again with trafo
N1.traf <- N1; trafo(N1.traf) <- tfct
N1R.traf <- N1.Rob; trafo(N1R.traf) <- tfct
IC1.traf <- optIC(model = N1R.traf, risk = asMSE())
(est0.traf <- MDEstimator(x, N1.traf))
(est1.traf <- kStepEstimator(x, IC1.traf, est0, steps = 3))
# or simply
(est2.traf <- oneStepEstimator(x, IC1.traf, est0))

### main: location; nuisance: scale
N1.NS <- L2LocationUnknownScaleFamily()
N1R.NS <- InfRobModel(center = N1.NS, neighbor = ContNeighborhood(radius = 0.5))
IC1.NS <- optIC(model = N1.NS, risk = asCov())
IC2.NS <- optIC(model = N1R.NS, risk = asMSE())
(est0.NS <- MDEstimator(x, N1.NS))
(est1.NS <- kStepEstimator(x, IC2.NS, est0, steps = 3))
(est2.NS <- oneStepEstimator(x, IC2.NS, est0))


## a simple example
library(MASS)
data(chem)
initial.est <- c(median(chem), mad(chem))
system.time(ROest1 <- roptest(chem, NormLocationScaleFamily(), eps.upper = 0.1, steps = 3L, 
                              initial.est = initial.est, useLast = TRUE))

## optimized for speed
library(RobLox)
system.time(ROest2 <- roblox(chem, eps.upper = 0.1, k = 3, returnIC = TRUE))

## comparison
estimate(ROest1)
estimate(ROest2)

## confidence intervals
confint(ROest1, symmetricBias())
confint(ROest2, symmetricBias())

########### again with trafo
system.time(ROest1.traf <- roptest(chem, N1.traf, eps.upper = 0.1, steps = 3L,
                              initial.est = initial.est, useLast = TRUE))
estimate(ROest1.traf)
confint(ROest1.traf, symmetricBias())
