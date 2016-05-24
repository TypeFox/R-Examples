###############################################################################
## Example: Gumbel Location Family
## computations numerically less stable than in case of the 
## Exponential Scale Family
###############################################################################
require(ROptEst)
options("newDevice"=TRUE)

## generates Gumbel Location Family with loc = 0
## (known scale = 1)
distrExOptions(ElowerTruncQuantile = 1e-7) # non-finite function value in integrate
G0 <- GumbelLocationFamily(loc = 0, scale = 1)
G0        # show G0
plot(G0)  # plot of Gumbel(loc = 0, scale = 1) and L_2 derivative
checkL2deriv(G0)

## classical optimal IC
G0.IC0 <- optIC(model = G0, risk = asCov())
G0.IC0       # show IC
plot(G0.IC0) # plot IC
checkIC(G0.IC0)

## L_2 family + infinitesimal neighborhood
G0.Rob1 <- InfRobModel(center = G0, neighbor = ContNeighborhood(radius = 0.5))
G0.Rob1     # show G0.Rob1
G0.Rob2 <- InfRobModel(center = G0, neighbor = TotalVarNeighborhood(radius = 0.5))

## MSE solution
E1.Rob1 <- InfRobModel(center = ExpScaleFamily(), neighbor = ContNeighborhood(radius = 0.5))
(E1.IC1 <- optIC(model=E1.Rob1, risk=asMSE()))
G0.IC1 <- optIC(model=G0.Rob1, risk=asMSE())
checkIC(G0.IC1)
Risks(G0.IC1)
clip(E1.IC1)
cent(E1.IC1)
stand(E1.IC1)
clip(G0.IC1)
cent(G0.IC1)
stand(G0.IC1)

## alternatively
G0.IC11 <- E1.IC1 # rate = 1!
CallL2Fam(G0.IC11) <- call("GumbelLocationFamily")
cent(G0.IC11) <- -cent(E1.IC1)
G0.IC11
checkIC(G0.IC11)
Risks(G0.IC11)

E1.Rob2 <- InfRobModel(center = ExpScaleFamily(), neighbor = TotalVarNeighborhood(radius = 0.5))
E1.IC2 <- optIC(model=E1.Rob2, risk=asMSE())
G0.IC2 <- optIC(model=G0.Rob2, risk=asMSE())
checkIC(G0.IC2)
Risks(G0.IC2)
clipLo(E1.IC2)
clipUp(E1.IC2)
stand(E1.IC2)
clipLo(G0.IC2)
clipUp(G0.IC2)
stand(G0.IC2)
## alternatively
G0.IC21 <- E1.IC2 # rate = 1!
CallL2Fam(G0.IC21) <- call("GumbelLocationFamily")
clipLo(G0.IC21) <- -clipUp(E1.IC2)
clipUp(G0.IC21) <- -clipLo(E1.IC2)
G0.IC21
checkIC(G0.IC21)
Risks(G0.IC21)

## lower case solutions
(G0.IC3 <- optIC(model=G0.Rob1, risk=asBias()))
checkIC(G0.IC3)
Risks(G0.IC3)
(G0.IC4 <- optIC(model=G0.Rob2, risk=asBias()))
checkIC(G0.IC4)
Risks(G0.IC4)

## Hampel solution
(G0.IC5 <- optIC(model=G0.Rob1, risk=asHampel(bound=clip(G0.IC1))))
checkIC(G0.IC5)
Risks(G0.IC5)
(G0.IC6 <- optIC(model=G0.Rob2, risk=asHampel(bound=Risks(G0.IC2)$asBias$value), maxiter = 100))
checkIC(G0.IC6)
Risks(G0.IC6)

## radius minimax IC
## numerically instable for small 'loRad'!
## => use connection to ExpScaleFamily for computations
(G0.IC7 <- radiusMinimaxIC(L2Fam=G0, neighbor=ContNeighborhood(), 
                risk=asMSE(), loRad=0.5, upRad=1.0))
checkIC(G0.IC7)
Risks(G0.IC7)
(G0.IC8 <- radiusMinimaxIC(L2Fam=G0, neighbor=TotalVarNeighborhood(), 
                risk=asMSE(), loRad=0.5, upRad=1.5))
checkIC(G0.IC8)
Risks(G0.IC8)

## least favorable radius
## numerically instable!
## => use connection to ExpScaleFamily for computations
(G0.r.rho1 <- leastFavorableRadius(L2Fam=G0, neighbor=ContNeighborhood(),
                    risk=asMSE(), rho=0.5))
(G0.r.rho2 <- leastFavorableRadius(L2Fam=G0, neighbor=TotalVarNeighborhood(),
                    risk=asMSE(), rho=1/3))

## one-step estimation
## 1. generate a contaminated sample
ind <- rbinom(1e2, size=1, prob=0.05) 
G0.x <- rgumbel(1e2, loc=(1-ind)*0.5+ind*1)

## 2. Kolmogorov(-Smirnov) minimum distance estimator
(G0.est01 <- MDEstimator(x=G0.x, GumbelLocationFamily()))

## 3. Cramer von Mises minimum distance estimator
(G0.est02 <- MDEstimator(x=G0.x, GumbelLocationFamily(), distance = CvMDist))

## 4. one-step estimation: radius known
G0.Rob3 <- InfRobModel(center=GumbelLocationFamily(loc=estimate(G0.est02)), 
                       neighbor=ContNeighborhood(radius=0.5))
G0.IC9 <- optIC(model=G0.Rob3, risk=asMSE())
(G0.est1 <- oneStepEstimator(G0.x, IC=G0.IC9, start=G0.est02))

## 5. k-step estimation: radius known
(G0.est2 <- kStepEstimator(G0.x, IC=G0.IC9, start=G0.est02, steps = 3))

## 6. M estimation: radius known
G0.Rob31 <- InfRobModel(center=GumbelLocationFamily(loc=0), 
                        neighbor=ContNeighborhood(radius=0.5))
G0.IC91 <- optIC(model=G0.Rob31, risk=asMSE())
(G0.est11 <- locMEstimator(G0.x, IC=G0.IC91))

## 7. one-step estimation: radius interval
G0.IC10 <- radiusMinimaxIC(L2Fam=GumbelLocationFamily(loc=estimate(G0.est02)),
                neighbor=ContNeighborhood(), risk=asMSE(), loRad=0.5, upRad=1)
(G0.est3 <- oneStepEstimator(G0.x, IC=G0.IC10, start=G0.est02))

## 8. k-step estimation: radius interval
(G0.est4 <- kStepEstimator(G0.x, IC=G0.IC10, start=G0.est02, steps = 3))

## 9. M estimation: radius interval
G0.IC101 <- radiusMinimaxIC(L2Fam=GumbelLocationFamily(),
                neighbor=ContNeighborhood(), risk=asMSE(), loRad=0.5, upRad=1)
(G0.est21 <- locMEstimator(G0.x, IC=G0.IC101))

## 10. It's easier to use function roptest for k-step (k >= 1) estimation!
sqrtn <- sqrt(length(G0.x))
G0.est5 <- roptest(G0.x, eps.lower = 0.5/sqrtn, eps.upper = 1/sqrtn, 
                   L2Fam = GumbelLocationFamily())
G0.est6 <- roptest(G0.x, eps.lower = 0.5/sqrtn, eps.upper = 1/sqrtn, 
                   L2Fam = GumbelLocationFamily(), steps = 3)

## comparison - radius known
estimate(G0.est1)
estimate(G0.est2)
estimate(G0.est11)

## confidence intervals
confint(G0.est1, symmetricBias())
confint(G0.est2, symmetricBias())
confint(G0.est11, symmetricBias())

## comparison - radius interval
estimate(G0.est3)
estimate(G0.est5)
estimate(G0.est4)
estimate(G0.est6)
estimate(G0.est21)

## confidence intervals
confint(G0.est3, symmetricBias())
confint(G0.est5, symmetricBias())
confint(G0.est4, symmetricBias())
confint(G0.est6, symmetricBias())
confint(G0.est21, symmetricBias())

distrExOptions(ElowerTruncQuantile=0) # default
