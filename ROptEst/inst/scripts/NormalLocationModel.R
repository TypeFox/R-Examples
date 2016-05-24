###############################################################################
## Example: Normal Location
###############################################################################
require(ROptEst)
options("newDevice"=TRUE)

## generates Normal Location Family with scale = 3
N0 <- NormLocationFamily(mean=2, sd=3)
N0        # show G0
plot(N0)  # plot of Norm(mean = 2, sd = 3) and L_2 derivative
checkL2deriv(N0)

## classical optimal IC
N0.IC0 <- optIC(model = N0, risk = asCov())
N0.IC0       # show IC
plot(N0.IC0) # plot IC
checkIC(N0.IC0)
Risks(N0.IC0)

## L_2 family + infinitesimal neighborhood
N0.Rob1 <- InfRobModel(center = N0, neighbor = ContNeighborhood(radius = 0.5))
N0.Rob1     # show N0.Rob1
N0.Rob2 <- InfRobModel(center = N0, neighbor = TotalVarNeighborhood(radius = 0.5))

## OBRE solution (ARE 95%)

system.time(N0.ICA <- optIC(model = N0.Rob1, risk = asAnscombe(.95),upper=NULL,lower=NULL, verbose=TRUE))
checkIC(N0.ICA)
Risks(N0.ICA)
plot(N0.ICA)

system.time(N0.ICA2 <- optIC(model = N0.Rob2, risk = asAnscombe(.95),upper=NULL,lower=NULL, verbose=TRUE))
checkIC(N0.ICA2)
Risks(N0.ICA2)
plot(N0.ICA2)


## MSE solution
(N0.IC1 <- optIC(model=N0.Rob1, risk=asMSE()))
checkIC(N0.IC1)
Risks(N0.IC1)
plot(N0.IC1)

(N0.IC2 <- optIC(model=N0.Rob2, risk=asMSE()))
checkIC(N0.IC2)
Risks(N0.IC2)
plot(N0.IC2)

## L1 solution

(N0.IC1.L1 <- optIC(model=N0.Rob1, risk=asL1()))
checkIC(N0.IC1.L1)
Risks(N0.IC1.L1)
plot(N0.IC1.L1)

(N0.IC2.L1 <- optIC(model=N0.Rob2, risk=asL1()))
checkIC(N0.IC2.L1)
Risks(N0.IC2.L1)
plot(N0.IC2.L1)

## L4 solution

(N0.IC1.L4 <- optIC(model=N0.Rob1, risk=asL4()))
checkIC(N0.IC1.L4)
Risks(N0.IC1.L4)
plot(N0.IC1.L4)

(N0.IC2.L4 <- optIC(model=N0.Rob2, risk=asL4()))
checkIC(N0.IC2.L4)
Risks(N0.IC2.L4)
plot(N0.IC2.L4)

## semivar-solution

(N0.IC1.SemiV.p <- optIC(model=N0.Rob1, risk=asSemivar(sign=1)))
checkIC(N0.IC1.SemiV.p)
Risks(N0.IC1.SemiV.p)
plot(N0.IC1.SemiV.p)

(N0.IC1.SemiV.m <- optIC(model=N0.Rob1, risk=asSemivar(sign=-1)))
checkIC(N0.IC1.SemiV.m)
Risks(N0.IC1.SemiV.m)
plot(N0.IC1.SemiV.m)

## lower case solutions
(N0.IC3 <- optIC(model=N0.Rob1, risk=asBias()))
checkIC(N0.IC3)
Risks(N0.IC3)
plot(N0.IC3)
(N0.IC4 <- optIC(model=N0.Rob2, risk=asBias()))
checkIC(N0.IC4)
Risks(N0.IC4)
plot(N0.IC4)

## Hampel solution
(N0.IC5 <- optIC(model=N0.Rob1, risk=asHampel(bound=clip(N0.IC1))))
checkIC(N0.IC5)
Risks(N0.IC5)
plot(N0.IC5)
(N0.IC6 <- optIC(model=N0.Rob2, risk=asHampel(bound=Risks(N0.IC2)$asBias$value), maxiter = 200))
checkIC(N0.IC6)
Risks(N0.IC6)
plot(N0.IC6)

## radius minimax IC
(N0.IC7 <- radiusMinimaxIC(L2Fam=N0, neighbor=ContNeighborhood(),
                risk=asMSE(), loRad=0, upRad=Inf))
checkIC(N0.IC7)
Risks(N0.IC7)
plot(N0.IC7)
(N0.IC8 <- radiusMinimaxIC(L2Fam=N0, neighbor=TotalVarNeighborhood(),
                risk=asMSE(), loRad=0, upRad=Inf))
checkIC(N0.IC8)
Risks(N0.IC8)
plot(N0.IC8)

neighbor.0 <- ContNeighborhood()
getReq(asMSE(),neighbor.0,N0.ICA,N0.IC1,n=1)
getReq(asMSE(),neighbor.0,N0.ICA,N0.IC1,n=30)
getReq(asL1(),neighbor.0,N0.ICA,N0.IC1,n=30)
getReq(asL4(),neighbor.0,N0.ICA,N0.IC1,n=30)
getReq(asMSE(),neighbor.0,N0.ICA,N0.IC7,n=30)
getReq(asL1(),neighbor.0,N0.ICA,N0.IC7,n=30)
getReq(asL4(),neighbor.0,N0.ICA,N0.IC7,n=30)
getReq(asMSE(),neighbor.0,N0.IC1,N0.IC7,n=30)


getMaxIneff(N0.ICA,neighbor=ContNeighborhood())
getMaxIneff(N0.IC1,neighbor=ContNeighborhood())
getMaxIneff(N0.IC7,neighbor=ContNeighborhood())



## least favorable radius
## (may take quite some time!)
(N0.r.rho1 <- leastFavorableRadius(L2Fam=N0, neighbor=ContNeighborhood(),
                    risk=asMSE(), rho=0.5))
(N0.r.rho2 <- leastFavorableRadius(L2Fam=N0, neighbor=TotalVarNeighborhood(),
                    risk=asMSE(), rho=1/3))


## For estimation use function roptest
ind <- rbinom(1e2, size=1, prob=0.05)
x <- rnorm(1e2, mean=1, sd = (1-ind)+ind*9)

## 1-step: contamination known
est1 <- roptest(x, eps = 0.05, L2Fam = NormLocationFamily(mean = 1))
est1v <- roptest(x, eps = 0.025, L2Fam = NormLocationFamily(mean = 1),
                 neighbor = TotalVarNeighborhood())

## k-step: contamination known
est2 <- roptest(x, eps = 0.05, L2Fam = NormLocationFamily(mean = 1), steps = 3)
est2v <- roptest(x, eps = 0.025, L2Fam = NormLocationFamily(mean = 1),
                 neighbor = TotalVarNeighborhood(), steps = 3)

## comparison
estimate(est1)
estimate(est2)
estimate(est1v)
estimate(est2v)

## confidence intervals
confint(est1, symmetricBias())
confint(est2, symmetricBias())
confint(est1v, symmetricBias())
confint(est2v, symmetricBias())
