###############################################################################
## Example: Normal Location
###############################################################################
require(ROptEst)
options("newDevice"=TRUE)

## generates Normal Location Family with mean = 0
N0 <- NormLocationFamily(mean=0) 
n <- 100
tau <- qnorm(0.975)

## classical optimal IC (radius = 0!)
N0.Rob1 <- InfRobModel(center = N0, neighbor = ContNeighborhood(radius = 0))
N0.Rob2 <- InfRobModel(center = N0, neighbor = TotalVarNeighborhood(radius = 0))

system.time(IC0c <- optIC(model=N0.Rob1, risk=asUnOvShoot(width = tau)))
checkIC(IC0c)
Risks(IC0c)
system.time(IC0v <- optIC(model=N0.Rob2, risk=asUnOvShoot(width = tau)))
checkIC(IC0v)
Risks(IC0v)

## boundary case
N0.Rob1 <- InfRobModel(center = N0, neighbor = ContNeighborhood(radius = 2*tau*1/sqrt(2*pi)))
N0.Rob2 <- InfRobModel(center = N0, neighbor = TotalVarNeighborhood(radius = tau*1/sqrt(2*pi)))

system.time(IC0c <- optIC(model=N0.Rob1, risk=asUnOvShoot(width = tau)))
checkIC(IC0c)
Risks(IC0c)
system.time(IC0v <- optIC(model=N0.Rob2, risk=asUnOvShoot(width = tau)))
checkIC(IC0v)
Risks(IC0v)


## L_2 family + infinitesimal resp. fixed neighborhood
rad <- 0.5
N0.Rob1 <- InfRobModel(center = N0, neighbor = ContNeighborhood(radius = rad))
N0.Rob2 <- InfRobModel(center = N0, neighbor = TotalVarNeighborhood(radius = rad/2))
N0.Rob3 <- FixRobModel(center = N0, neighbor = ContNeighborhood(radius = rad/sqrt(n)))
N0.Rob4 <- FixRobModel(center = N0, neighbor = TotalVarNeighborhood(radius = rad/2/sqrt(n)))

## asUnOvShoot solution
N0.IC1 <- optIC(model = N0.Rob1, risk = asUnOvShoot(width = tau))
checkIC(N0.IC1)
Risks(N0.IC1)
plot(N0.IC1)

N0.IC2 <- optIC(model = N0.Rob2, risk = asUnOvShoot(width = tau))
checkIC(N0.IC2)
Risks(N0.IC2)
plot(N0.IC2)

## fiUnOvShoot solution
N0.IC3 <- optIC(model=N0.Rob3, risk=fiUnOvShoot(width = tau/sqrt(n)), sampleSize = n)
checkIC(N0.IC3)
Risks(N0.IC3)
plot(N0.IC3)

N0.IC4 <- optIC(model=N0.Rob4, risk=fiUnOvShoot(width = tau/sqrt(n)), sampleSize = n)
checkIC(N0.IC4)
Risks(N0.IC4)
plot(N0.IC4)

## O(n^(-0.5))-corrected solution 
## in case of contamination neighborhoods
N0.IC5 <- N0.IC1
clipUp1 <- clipUp(N0.IC1)/as.vector(stand(N0.IC1))
clipUp5 <- max(0, clipUp1 - rad*(rad + clipUp1*tau)/(sqrt(n)*2*tau*pnorm(-clipUp1)))
stand5 <- 1/(2*pnorm(clipUp5)-1)
clipUp(N0.IC5) <- stand5*clipUp5
clipLo(N0.IC5) <- -clipUp(N0.IC5)
stand(N0.IC5) <- as.matrix(stand5)
Infos(N0.IC5) <- matrix(c("manually", "O(n^(-1/2))-corrected solution"), ncol=2,
                    dimnames=list(character(0), c("method", "message")))
checkIC(N0.IC5)
getRiskIC(N0.IC5, asUnOvShoot(width = tau), ContNeighborhood(radius=rad))
getRiskIC(N0.IC5, fiUnOvShoot(width = tau/sqrt(n)), ContNeighborhood(radius=rad/sqrt(n)), sampleSize = n)

## O(n^(-1))-corrected solution 
## in case of total variation neighborhoods
N0.IC6 <- N0.IC2
clipUp2 <- clipUp(N0.IC2)/as.vector(stand(N0.IC2))
clipUp6 <- max(0, clipUp2 - tau*(2*clipUp2^2*rad/2 + tau*dnorm(clipUp2))/(6*n*pnorm(-clipUp2)))
stand6 <- 1/(2*pnorm(clipUp6)-1)
clipUp(N0.IC6) <- stand6*clipUp6
clipLo(N0.IC6) <- -clipUp(N0.IC6)
stand(N0.IC6) <- as.matrix(stand6)
Infos(N0.IC6) <- matrix(c("manually", "O(n^(-1))-corrected solution"), ncol=2,
                    dimnames=list(character(0), c("method", "message")))
checkIC(N0.IC6)
getRiskIC(N0.IC6, asUnOvShoot(width = tau), TotalVarNeighborhood(radius=rad/2))
getRiskIC(N0.IC6, fiUnOvShoot(width = tau/sqrt(n)), TotalVarNeighborhood(radius=rad/2/sqrt(n)), sampleSize = n)


## estimation
## 1. generate a contaminated sample
ind <- rbinom(1e2, size = 1, prob = 0.05) 
X <- rnorm(1e2, mean = ind*4)
summary(X)

## 2. M estimation
N0.Rob5 <- InfRobModel(center = NormLocationFamily(mean = 0), 
                       neighbor = ContNeighborhood(radius = 0.5))
N0.IC7 <- optIC(model=N0.Rob5, risk=asUnOvShoot(width = 1.960))
(Mest1 <- locMEstimator(X, IC=N0.IC7))

## confidence interval
confint(Mest1, symmetricBias())

N0.Rob6 <- FixRobModel(center = NormLocationFamily(mean = 0), 
                       neighbor = ContNeighborhood(radius = 0.05))
N0.IC8 <- optIC(model = N0.Rob6, risk=fiUnOvShoot(width = 1.960/sqrt(n)), sampleSize = 1e2)
(Mest2 <- locMEstimator(X, IC=N0.IC8))

## confidence interval
confint(Mest2, symmetricBias())

## 3. Kolmogorov(-Smirnov) minimum distance estimator
(est0 <- MDEstimator(x=X, NormLocationFamily()))

## 4. one-step estimation
N0.Rob7 <- InfRobModel(center = NormLocationFamily(mean = estimate(est0)), 
                       neighbor = ContNeighborhood(radius=0.5))
N0.IC9 <- optIC(model=N0.Rob7, risk=asUnOvShoot(width = 1.960))
(est1 <- oneStepEstimator(X, IC = N0.IC9, start = est0))
confint(est1, symmetricBias())

N0.Rob8 <- FixRobModel(center = NormLocationFamily(mean = estimate(est0)), 
                       neighbor = ContNeighborhood(radius=0.05))
N0.IC10 <- optIC(model=N0.Rob8, risk=fiUnOvShoot(width = 1.960/sqrt(n)), sampleSize = 1e2)
(est2 <- oneStepEstimator(X, IC = N0.IC10, start = est0))
confint(est2, symmetricBias())
