###############################################################################
## Example: Poisson Family
###############################################################################
require(ROptEst)
options("newDevice"=TRUE)

distroptions("TruncQuantile"= 1e-10) # increases numerical support of Pois;
                                     # i.e., increases precision of the 
                                     # computations
## generates Poisson Family with theta = 4.5
P <- PoisFamily(lambda = 4.5) 
P       # show P
plot(P) # plot of Pois(lambda = 4.5) and L_2 derivative
checkL2deriv(P)

## classical optimal IC
IC0 <- optIC(model = P, risk = asCov())
IC0       # show IC
checkIC(IC0)
Risks(IC0)
plot(IC0) # plot IC

## L_2 family + infinitesimal neighborhood
RobP1 <- InfRobModel(center = P, neighbor = ContNeighborhood(radius = 0.5))
RobP1     # show RobP1
(RobP2 <- InfRobModel(center = P, neighbor = TotalVarNeighborhood(radius = 0.5)))

## lower case radius
lowerCaseRadius(L2Fam = P, ContNeighborhood(radius = 0.5), risk = asMSE())
lowerCaseRadius(L2Fam = P, TotalVarNeighborhood(radius = 0.5), risk = asMSE())

## OBRE solution (ARE = .95)
system.time(ICA <- optIC(model = RobP1, risk = asAnscombe(.95),
                         upper=NULL,lower=NULL, verbose=TRUE))
checkIC(ICA)
Risks(ICA)
plot(ICA)

system.time(ICA.p <- optIC(model = RobP1, 
                         risk = asAnscombe(.95,biastype=positiveBias()),
                         upper=NULL,lower=NULL, verbose=TRUE))
checkIC(ICA.p)
Risks(ICA.p)
plot(ICA.p)

## MSE solution
(IC1 <- optIC(model=RobP1, risk=asMSE()))
checkIC(IC1)
Risks(IC1)
plot(IC1)

(IC1.p <- optIC(model=RobP1, risk=asMSE(biastype=positiveBias())))
checkIC(IC1.p)
Risks(IC1.p)
plot(IC1.p)

(IC1.a <- optIC(model=RobP1, risk=asMSE(biastype=asymmetricBias(nu = c(1,0.2)))))
checkIC(IC1.a)
Risks(IC1.a)
plot(IC1.a)

(IC2 <- optIC(model=RobP2, risk=asMSE()))
checkIC(IC2)
Risks(IC2)
plot(IC2)


## lower case solutions
(IC3 <- optIC(model=RobP1, risk=asBias()))
checkIC(IC3)
Risks(IC3)
plot(IC3)

(IC3.p <- optIC(model=RobP1, risk=asBias(biastype=positiveBias())))
checkIC(IC3.p) # numerical problem???
Risks(IC3.p)
plot(IC3.p)

(IC3.a <- optIC(model=RobP1, risk=asBias(biastype=asymmetricBias(nu = c(1,0.2)))))
checkIC(IC3.a)
Risks(IC3.a)
plot(IC3.a)


(IC4 <- optIC(model=RobP2, risk=asBias()))
checkIC(IC4)
Risks(IC4)
plot(IC4)

## Hampel solution
(IC5 <- optIC(model=RobP1, risk=asHampel(bound=clip(IC1))))
checkIC(IC5)
Risks(IC5)
plot(IC5)

(IC6 <- optIC(model=RobP2, risk=asHampel(bound=Risks(IC2)$asBias$value), maxiter = 200))
checkIC(IC6)
Risks(IC6)
plot(IC6)


## radius minimax IC
(IC7 <- radiusMinimaxIC(L2Fam=P, neighbor=ContNeighborhood(), 
                risk=asMSE(), loRad=0, upRad=0.5))
checkIC(IC7)
Risks(IC7)
plot(IC7)

(IC7.p <- radiusMinimaxIC(L2Fam=P, neighbor=ContNeighborhood(), 
                risk=asMSE(biastype=positiveBias()), loRad=0, upRad=0.5))
checkIC(IC7.p)
Risks(IC7.p)
plot(IC7.p)

(IC7.a <- radiusMinimaxIC(L2Fam=P, neighbor=ContNeighborhood(), 
                risk=asMSE(biastype=asymmetricBias(nu = c(1,0.2))), loRad=0, upRad=0.5))
checkIC(IC7.a)
Risks(IC7.a)
plot(IC7.a)

(IC8 <- radiusMinimaxIC(L2Fam=P, neighbor=TotalVarNeighborhood(), 
                risk=asMSE(), loRad=0, upRad=Inf))
checkIC(IC8)
Risks(IC8)
plot(IC8)

## least favorable radius
## (may take quite some time!)
(r.rho1 <- leastFavorableRadius(L2Fam=P, neighbor=ContNeighborhood(),
                    risk=asMSE(), rho=0.5))
(r.rho2 <- leastFavorableRadius(L2Fam=P, neighbor=TotalVarNeighborhood(),
                    risk=asMSE(), rho=1/3))

## one-step estimation
## Example: Rutherford-Geiger (1910)
## cf. Feller~(1968), Section VI.7 (a)
x <- c(rep(0, 57), rep(1, 203), rep(2, 383), rep(3, 525), rep(4, 532), 
       rep(5, 408), rep(6, 273), rep(7, 139), rep(8, 45), rep(9, 27), 
       rep(10, 10), rep(11, 4), rep(12, 0), rep(13, 1), rep(14, 1))

## 0. ML-estimator (mean)
(est0 <- mean(x))

## with distrMod
MLEstimator(x, PoisFamily())

## with MASS
library(MASS)
fitdistr(x, densfun = "Poisson")

## 1.1. Kolmogorov(-Smirnov) minimum distance estimator
(est11 <- MDEstimator(x=x, PoisFamily()))

## 1.2. Cramer von Mises minimum distance estimator
(est12 <- MDEstimator(x=x, PoisFamily(), distance = CvMDist))

## 2. k-step estimation: contamination neighborhood
## 2.1 small amount of contamination < 2%
IC9 <- radiusMinimaxIC(L2Fam=PoisFamily(lambda=estimate(est11)),
                neighbor=ContNeighborhood(), risk=asMSE(), loRad=0, upRad=1)
IC10 <- radiusMinimaxIC(L2Fam=PoisFamily(lambda=estimate(est12)),
                neighbor=ContNeighborhood(), risk=asMSE(), loRad=0, upRad=1)
(est211 <- kStepEstimator(x, IC=IC9, start=est11, steps = 1L))
## one could also use function oneStepEstimator
oneStepEstimator(x, IC=IC9, start=est11)
checkIC(pIC(est211))

(est212 <- kStepEstimator(x, IC=IC9, start=est11, steps = 3L))
checkIC(pIC(est212))

(est213 <- kStepEstimator(x, IC=IC10, start=est12, steps = 1L))
checkIC(pIC(est213))

(est214 <- kStepEstimator(x, IC=IC10, start=est12, steps = 3L))
checkIC(pIC(est214))

## Its simpler to use roptest!
(est215 <- roptest(x, PoisFamily(), eps.upper = 1/sqrt(length(x)), steps = 3L))
checkIC(pIC(est215))

## comparision of estimates
estimate(est211)
estimate(est212)
estimate(est213)
estimate(est214)
estimate(est215)

## confidence intervals
confint(est211, symmetricBias())
confint(est212, symmetricBias())
confint(est213, symmetricBias())
confint(est214, symmetricBias())
confint(est215, symmetricBias())


## 2.2 amount of contamination unknown
IC11 <- radiusMinimaxIC(L2Fam=PoisFamily(lambda=estimate(est11)),
                neighbor=ContNeighborhood(), risk=asMSE(), loRad=0, upRad=Inf)
IC12 <- radiusMinimaxIC(L2Fam=PoisFamily(lambda=estimate(est12)),
                neighbor=ContNeighborhood(), risk=asMSE(), loRad=0, upRad=Inf)
(est221 <- oneStepEstimator(x, IC=IC11, start=est11))
kStepEstimator(x, IC=IC11, start=est11)
checkIC(pIC(est221))

(est222 <- kStepEstimator(x, IC=IC11, start=est11, steps = 3L))
checkIC(pIC(est222))

(est223 <- kStepEstimator(x, IC=IC12, start=est12, steps = 1L))
checkIC(pIC(est223))

(est224 <- kStepEstimator(x, IC=IC12, start=est12, steps = 3L))
checkIC(pIC(est224))

## Its simpler to use roptest!
(est225 <- roptest(x, PoisFamily(), eps.upper = 0.5, steps = 3L))
checkIC(pIC(est225))

## comparision of estimates
estimate(est221)
estimate(est222)
estimate(est223)
estimate(est224)
estimate(est225)

## confidence intervals
confint(est221, symmetricBias())
confint(est222, symmetricBias())
confint(est223, symmetricBias())
confint(est224, symmetricBias())
confint(est225, symmetricBias())

## 3. k-step estimation: total variation neighborhood
## 3.1 small amount of contamination < 2%
IC13 <- radiusMinimaxIC(L2Fam=PoisFamily(lambda=estimate(est11)),
                neighbor=TotalVarNeighborhood(), risk=asMSE(), loRad=0, upRad=1)
IC14 <- radiusMinimaxIC(L2Fam=PoisFamily(lambda=estimate(est12)),
                neighbor=TotalVarNeighborhood(), risk=asMSE(), loRad=0, upRad=1)
(est311 <- kStepEstimator(x, IC=IC13, start=est11, steps = 1L))
## one could also use function oneStepEstimator
oneStepEstimator(x, IC=IC13, start=est11)
checkIC(pIC(est311))

(est312 <- kStepEstimator(x, IC=IC13, start=est11, steps = 3L))
checkIC(pIC(est312))

(est313 <- kStepEstimator(x, IC=IC14, start=est12, steps = 1L))
checkIC(pIC(est313))

(est314 <- kStepEstimator(x, IC=IC14, start=est12, steps = 3L))
checkIC(pIC(est314))

## Its simpler to use roptest!
(est315 <- roptest(x, PoisFamily(), eps.upper = 1/sqrt(length(x)), steps = 3L, 
                   neighbor = TotalVarNeighborhood()))
checkIC(pIC(est315))

(est316 <- kStepEstimator(x, IC=IC14, start=est12, steps = 15L))
ksteps(e316)
ksteps(e316, diff = TRUE)


## comparison of estimates
estimate(est311)
estimate(est312)
estimate(est313)
estimate(est314)
estimate(est315)

## confidence intervals
confint(est311, symmetricBias())
confint(est312, symmetricBias())
confint(est313, symmetricBias())
confint(est314, symmetricBias())
confint(est315, symmetricBias())


## 3.2 amount of contamination unknown
IC15 <- radiusMinimaxIC(L2Fam=PoisFamily(lambda=estimate(est11)),
                neighbor=TotalVarNeighborhood(), risk=asMSE(), loRad=0,
                upRad=Inf, loRad0 = 0.02, verbose = TRUE)
IC16 <- radiusMinimaxIC(L2Fam=PoisFamily(lambda=estimate(est12)),
                neighbor=TotalVarNeighborhood(), risk=asMSE(), loRad=0,
                upRad=Inf, loRad0 = 0.02, verbose = TRUE)
(est321 <- oneStepEstimator(x, IC=IC15, start=est11))
kStepEstimator(x, IC=IC15, start=est11)
checkIC(pIC(est321))

(est322 <- kStepEstimator(x, IC=IC15, start=est11, steps = 3L))
checkIC(pIC(est322))

(est323 <- kStepEstimator(x, IC=IC16, start=est12, steps = 1L))
checkIC(pIC(est323))

(est324 <- kStepEstimator(x, IC=IC16, start=est12, steps = 3L))
checkIC(pIC(est324))

## Its simpler to use roptest!
(est325 <- roptest(x, PoisFamily(), eps.upper = 0.5, steps = 3L, 
                   neighbor = TotalVarNeighborhood()))
checkIC(pIC(est325))

## comparision of estimates
estimate(est321)
estimate(est322)
estimate(est323)
estimate(est324)
estimate(est325)

## confidence intervals
confint(est321, symmetricBias())
confint(est322, symmetricBias())
confint(est323, symmetricBias())
confint(est324, symmetricBias())
confint(est325, symmetricBias())

distroptions("TruncQuantile"= 1e-5) # default
