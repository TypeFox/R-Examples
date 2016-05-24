###############################################################################
## Example: Lognormal Scale and Normal Location
###############################################################################
require(ROptEst)
options("newDevice"=TRUE)

## generates Lognormal Scale Family with meanlog = 0, sdlog = 1
LN1 <- LnormScaleFamily() 
LN1        # show LN1
plot(LN1)  # plot of Lnorm() and L_2 derivative
checkL2deriv(LN1)

## generates Normal Location Family with mean = 0
N0 <- NormLocationFamily(mean=0, sd=1) 
N0        # show G0
plot(N0)  # plot of Norm(mean = 0, sd = 1) and L_2 derivative
checkL2deriv(N0)


## classical optimal IC
LN1.IC0 <- optIC(model = LN1, risk = asCov())
LN1.IC0       # show IC
plot(LN1.IC0) # plot IC
checkIC(LN1.IC0)
Risks(LN1.IC0)
N0.IC0 <- optIC(model = N0, risk = asCov())
N0.IC0       # show IC
plot(N0.IC0) # plot IC
checkIC(N0.IC0)
Risks(N0.IC0)


## L_2 family + infinitesimal neighborhood
LN1.Rob1 <- InfRobModel(center = LN1, neighbor = ContNeighborhood(radius = 0.5))
LN1.Rob1     # show LN1.Rob1
LN1.Rob2 <- InfRobModel(center = LN1, neighbor = TotalVarNeighborhood(radius = 0.25))
N0.Rob1 <- InfRobModel(center = N0, neighbor = ContNeighborhood(radius = 0.5))
N0.Rob1     # show N0.Rob1
N0.Rob2 <- InfRobModel(center = N0, neighbor = TotalVarNeighborhood(radius = 0.25))


## MSE solution
LN1.IC1 <- optIC(model=LN1.Rob1, risk=asMSE())
checkIC(LN1.IC1)
Risks(LN1.IC1)
plot(LN1.IC1)

N0.IC1 <- optIC(model=N0.Rob1, risk=asMSE())
checkIC(N0.IC1)
Risks(N0.IC1)
plot(N0.IC1)

clip(LN1.IC1)
cent(LN1.IC1)
stand(LN1.IC1)
clip(N0.IC1)
cent(N0.IC1)
stand(N0.IC1)

LN1.IC2 <- optIC(model=LN1.Rob2, risk=asMSE())
checkIC(LN1.IC2)
Risks(LN1.IC2)
plot(LN1.IC2)

N0.IC2 <- optIC(model=N0.Rob2, risk=asMSE())
checkIC(N0.IC2)
Risks(N0.IC2)
plot(N0.IC2)

clipLo(LN1.IC2)
clipUp(LN1.IC2)
stand(LN1.IC2)
clipLo(N0.IC2)
clipUp(N0.IC2)
stand(N0.IC2)


## lower case solutions
LN1.IC3 <- optIC(model=LN1.Rob1, risk=asBias())
checkIC(LN1.IC3)
Risks(LN1.IC3)
plot(LN1.IC3)

N0.IC3 <- optIC(model=N0.Rob1, risk=asBias())
checkIC(N0.IC3)
Risks(N0.IC3)
plot(N0.IC3)

LN1.IC4 <- optIC(model=LN1.Rob2, risk=asBias())
checkIC(LN1.IC4)
Risks(LN1.IC4)
plot(LN1.IC4)

N0.IC4 <- optIC(model=N0.Rob2, risk=asBias())
checkIC(N0.IC4)
Risks(N0.IC4)
plot(N0.IC4)


## Hampel solution
LN1.IC5 <- optIC(model=LN1.Rob1, risk=asHampel(bound=clip(LN1.IC1)))
checkIC(LN1.IC5)
Risks(LN1.IC5)
plot(LN1.IC5)

N0.IC5 <- optIC(model=N0.Rob1, risk=asHampel(bound=clip(N0.IC1)))
checkIC(N0.IC5)
Risks(N0.IC5)
plot(N0.IC5)

LN1.IC6 <- optIC(model=LN1.Rob2, risk=asHampel(bound=Risks(LN1.IC2)$asBias$value))
checkIC(LN1.IC6)
Risks(LN1.IC6)
plot(LN1.IC6)

N0.IC6 <- optIC(model=N0.Rob2, risk=asHampel(bound=Risks(N0.IC2)$asBias$value))
checkIC(N0.IC6)
Risks(N0.IC6)
plot(N0.IC6)

## radius minimax IC
(LN1.IC7 <- radiusMinimaxIC(L2Fam=LN1, neighbor=ContNeighborhood(), 
                risk=asMSE(), loRad=0, upRad=0.5, loRad0=2e-3))
checkIC(LN1.IC7)
Risks(LN1.IC7)
plot(LN1.IC7)

(N0.IC7 <- radiusMinimaxIC(L2Fam=N0, neighbor=ContNeighborhood(), 
                risk=asMSE(), loRad=0, upRad=0.5, loRad0=2e-3))
checkIC(N0.IC7)
Risks(N0.IC7)
plot(N0.IC7)

(LN1.IC8 <- radiusMinimaxIC(L2Fam=LN1, neighbor=TotalVarNeighborhood(), 
                risk=asMSE(), loRad=0, upRad=0.25, loRad0=0.04))
checkIC(LN1.IC8)
Risks(LN1.IC8)
plot(LN1.IC8)

(N0.IC8 <- radiusMinimaxIC(L2Fam=N0, neighbor=TotalVarNeighborhood(), 
                risk=asMSE(), loRad=0, upRad=0.25, loRad0=0.04))
checkIC(N0.IC8)
Risks(N0.IC8)
plot(N0.IC8)


## least favorable radius
(LN1.r.rho1 <- leastFavorableRadius(L2Fam=LN1, neighbor=ContNeighborhood(),
                    risk=asMSE(), rho=0.5))
(N0.r.rho1 <- leastFavorableRadius(L2Fam=N0, neighbor=ContNeighborhood(),
                    risk=asMSE(), rho=0.5))
(LN1.r.rho2 <- leastFavorableRadius(L2Fam=LN1, neighbor=TotalVarNeighborhood(),
                    risk=asMSE(), rho=1/3))
(N0.r.rho2 <- leastFavorableRadius(L2Fam=N0, neighbor=TotalVarNeighborhood(),
                    risk=asMSE(), rho=1/3))

## For estimation use function roptest
ind <- rbinom(1e2, size=1, prob=0.05) 
x <- rnorm(1e2, mean=(1-ind)+ind*9)
y <- exp(x)

## 1-step: contamination known
est1 <- roptest(x, eps = 0.05, L2Fam = NormLocationFamily())
est2 <- roptest(y, eps = 0.05, L2Fam = LnormScaleFamily())

## k-step: contamination known
est3 <- roptest(x, eps = 0.05, L2Fam = NormLocationFamily(), steps = 3)
est4 <- roptest(y, eps = 0.05, L2Fam = LnormScaleFamily(), steps = 3)

## comparison
estimate(est1)
log(estimate(est2))
estimate(est3)
log(estimate(est4))

## confidence intervals
confint(est1, symmetricBias())
confint(est2, symmetricBias())
confint(est3, symmetricBias())
confint(est4, symmetricBias())
