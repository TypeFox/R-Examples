###############################################################################
# Normal linear regression
###############################################################################
system.time(require(ROptRegTS))

tau <- qnorm(0.995)
n <- 100

# Regressor distributions
K1 <- DiscreteDistribution(c(0.5, 1.0, 1.5))
K2 <- Unif(Min = -1, Max = 2)

# generate a normal linear regression family
LM1 <- NormLinRegFamily(RegDistr = K1)
LM2 <- NormLinRegFamily(RegDistr = K2)

## classical optimal IC (radius = 0!)
RobLM1c <- InfRobRegTypeModel(center = LM1, neighbor = ContNeighborhood(radius = 0))
RobLM1v <- InfRobRegTypeModel(center = LM1, neighbor = TotalVarNeighborhood(radius = 0))
RobLM2c <- InfRobRegTypeModel(center = LM2, neighbor = ContNeighborhood(radius = 0))
RobLM2v <- InfRobRegTypeModel(center = LM2, neighbor = TotalVarNeighborhood(radius = 0))

system.time(IC10c <- optIC(model=RobLM1c, risk=asUnOvShoot(width = tau)), gcFirst = TRUE)
checkIC(IC10c)
Risks(IC10c)
system.time(IC10v <- optIC(model=RobLM1v, risk=asUnOvShoot(width = tau)), gcFirst = TRUE)
checkIC(IC10v)
Risks(IC10v)
system.time(IC20c <- optIC(model=RobLM2c, risk=asUnOvShoot(width = tau)), gcFirst = TRUE)
checkIC(IC20c)
Risks(IC20c)
system.time(IC20v <- optIC(model=RobLM2v, risk=asUnOvShoot(width = tau)), gcFirst = TRUE)
checkIC(IC20v)
Risks(IC20v)

## boundary case
bound <- tau*E(K1, abs)/sqrt(2*pi)
RobLM1c <- InfRobRegTypeModel(center = LM1, neighbor = ContNeighborhood(radius = 2*bound))
RobLM1v <- InfRobRegTypeModel(center = LM1, neighbor = TotalVarNeighborhood(radius = bound))
bound <- tau*E(K2, abs)/sqrt(2*pi)
RobLM2c <- InfRobRegTypeModel(center = LM2, neighbor = ContNeighborhood(radius = 2*bound))
RobLM2v <- InfRobRegTypeModel(center = LM2, neighbor = TotalVarNeighborhood(radius = bound))

system.time(IC11c <- optIC(model=RobLM1c, risk=asUnOvShoot(width = tau)), gcFirst = TRUE)
checkIC(IC11c)
Risks(IC11c)
system.time(IC11v <- optIC(model=RobLM1v, risk=asUnOvShoot(width = tau)), gcFirst = TRUE)
checkIC(IC11v)
Risks(IC11v)
system.time(IC21c <- optIC(model=RobLM2c, risk=asUnOvShoot(width = tau)), gcFirst = TRUE)
checkIC(IC21c)
Risks(IC21c)
system.time(IC21v <- optIC(model=RobLM2v, risk=asUnOvShoot(width = tau)), gcFirst = TRUE)
checkIC(IC21v)
Risks(IC21v)


###############################################################################
## asUnOvShoot solution
###############################################################################
## infinitesimal robust model
rad <- 0.2
RobLM1c <- InfRobRegTypeModel(center = LM1, neighbor = ContNeighborhood(radius = rad))
RobLM1v <- InfRobRegTypeModel(center = LM1, neighbor = TotalVarNeighborhood(radius = rad/2))
RobLM2c <- InfRobRegTypeModel(center = LM2, neighbor = ContNeighborhood(radius = rad))
RobLM2v <- InfRobRegTypeModel(center = LM2, neighbor = TotalVarNeighborhood(radius = rad/2))

## asymptotic optimal IC
system.time(IC12c <- optIC(model=RobLM1c, risk=asUnOvShoot(width = tau)), gcFirst = TRUE)
checkIC(IC12c) # takes some time
Risks(IC12c)

system.time(IC12v <- optIC(model=RobLM1v, risk=asUnOvShoot(width = tau)), gcFirst = TRUE)
checkIC(IC12v) # takes some time
Risks(IC12v)

system.time(IC22c <- optIC(model=RobLM2c, risk=asUnOvShoot(width = tau)), gcFirst = TRUE)
checkIC(IC22c) # takes some time
Risks(IC22c)

system.time(IC22v <- optIC(model=RobLM2v, risk=asUnOvShoot(width = tau)), gcFirst = TRUE)
checkIC(IC22v) # takes some time
Risks(IC22v)

###############################################################################
## fiUnOvShoot solution
###############################################################################
## fixed robust model
RobLM3c <- FixRobRegTypeModel(center = LM1, neighbor = ContNeighborhood(radius = rad/sqrt(n)))
RobLM3v <- FixRobRegTypeModel(center = LM1, neighbor = TotalVarNeighborhood(radius = rad/2/sqrt(n)))
RobLM4c <- FixRobRegTypeModel(center = LM2, neighbor = ContNeighborhood(radius = rad/sqrt(n)))
RobLM4v <- FixRobRegTypeModel(center = LM2, neighbor = TotalVarNeighborhood(radius = rad/2/sqrt(n)))

system.time(IC32c <- optIC(model=RobLM3c, risk=fiUnOvShoot(width = tau/sqrt(n)), sampleSize = n), gcFirst = TRUE)
checkIC(IC32c) # takes some time
Risks(IC32c)

system.time(IC32v <- optIC(model=RobLM3v, risk=fiUnOvShoot(width = tau/sqrt(n)), sampleSize = n), gcFirst = TRUE)
checkIC(IC32v) # takes some time
Risks(IC32v)

system.time(IC42c <- optIC(model=RobLM4c, risk=fiUnOvShoot(width = tau/sqrt(n)), sampleSize = n), gcFirst = TRUE)
checkIC(IC42c) # takes some time
Risks(IC42c)

system.time(IC42v <- optIC(model=RobLM4v, risk=fiUnOvShoot(width = tau/sqrt(n)), sampleSize = n), gcFirst = TRUE)
checkIC(IC42v) # takes some time
Risks(IC42v)


# function for computation of standardizing constant
Afct <- function(x, c0){
    if(x == 0) return(0)
    return(x^2*pnorm(c0/abs(x)))
}

# O(n^(-0.5))-corrected solution 
# in case of contamination neighborhoods
IC13c <- IC12c
clipUp1 <- clipUp(IC12c)/as.vector(stand(IC12c))
clipUp3 <- max(0, clipUp1 - rad*(rad + clipUp1*tau)/(sqrt(n)*2*tau*E(K1, function(x, b){pnorm(-b/abs(x))}, b = clipUp1)))
stand3 <- 1/(2*E(K1, Afct, c0 = clipUp3) - E(K1, function(x){x^2}))
clipUp(IC13c) <- stand3*clipUp3
clipLo(IC13c) <- -clipUp(IC13c)
stand(IC13c) <- as.matrix(stand3)
checkIC(IC13c)
Risks(IC13c) <- getFiRiskRegTS(risk = fiUnOvShoot(width = tau/sqrt(n)), 
                   ErrorDistr = Norm(), Regressor = K1, 
                   neighbor = ContNeighborhood(radius = rad/sqrt(n)), 
                   clip = clipUp3, stand = stand3, sampleSize = n,
                   Algo = "A", cont = "left")
Infos(IC13c) <- matrix(c("manually", "O(n^(-1/2))-corrected solution"), ncol=2,
                    dimnames=list(character(0), c("method", "message")))
Risks(IC13c)
Risks(IC32c)

IC23c <- IC22c
clipUp1 <- clipUp(IC22c)/as.vector(stand(IC22c))
clipUp3 <- max(0, clipUp1 - rad*(rad + clipUp1*tau)/(sqrt(n)*2*tau*E(K1, function(x, b){pnorm(-b/abs(x))}, b = clipUp1)))
stand3 <- 1/(2*E(K2, Afct, c0 = clipUp3) - E(K2, function(x){x^2}))
clipUp(IC23c) <- stand3*clipUp3
clipLo(IC23c) <- -clipUp(IC23c)
stand(IC23c) <- as.matrix(stand3)
checkIC(IC23c)
Risks(IC23c) <- getFiRiskRegTS(risk = fiUnOvShoot(width = tau/sqrt(n)), 
                   ErrorDistr = Norm(), Regressor = K2, 
                   neighbor = ContNeighborhood(radius = rad/sqrt(n)), 
                   clip = clipUp3, stand = stand3, sampleSize = n,
                   Algo = "A", cont = "left")
Infos(IC23c) <- matrix(c("manually", "O(n^(-1/2))-corrected solution"), ncol=2,
                    dimnames=list(character(0), c("method", "message")))
Risks(IC23c)
Risks(IC42c)

# O(n^(-1))-corrected solution 
# in case of total variation neighborhoods
IC13v <- IC12v
clipUp1 <- clipUp(IC12v)/as.vector(stand(IC12v))
h1 <- E(K1, function(x, b){abs(x)^3*dnorm(b/abs(x))}, b = clipUp1)
h2 <- E(K1, function(x, b){pnorm(-b/abs(x))}, b = clipUp1)
clipUp3 <- max(0, clipUp1 - tau*(2*clipUp1^2*rad/2 + tau*h1)/(6*n*h2))
stand3 <- 1/(2*E(K1, Afct, c0 = clipUp3) - E(K1, function(x){x^2}))
clipUp(IC13v) <- stand3*clipUp3
clipLo(IC13v) <- -clipUp(IC13v)
stand(IC13v) <- as.matrix(stand3)
checkIC(IC13v)
Risks(IC13v) <- getFiRiskRegTS(risk = fiUnOvShoot(width = tau/sqrt(n)), 
                   ErrorDistr = Norm(), Regressor = K1, 
                   neighbor = TotalVarNeighborhood(radius = rad/2/sqrt(n)), 
                   clip = clipUp3, stand = stand3, sampleSize = n,
                   Algo = "A", cont = "left")
Infos(IC13v) <- matrix(c("manually", "O(n^(-1))-corrected solution"), ncol=2,
                    dimnames=list(character(0), c("method", "message")))
Risks(IC13v)
Risks(IC32v)

IC23v <- IC22v
clipUp1 <- clipUp(IC22v)/as.vector(stand(IC22v))
h1 <- E(K1, function(x, b){abs(x)^3*dnorm(b/abs(x))}, b = clipUp1)
h2 <- E(K1, function(x, b){pnorm(-b/abs(x))}, b = clipUp1)
clipUp3 <- max(0, clipUp1 - tau*(2*clipUp1^2*rad/2 + tau*h1)/(6*n*h2))
stand3 <- 1/(2*E(K2, Afct, c0 = clipUp3) - E(K2, function(x){x^2}))
clipUp(IC23v) <- stand3*clipUp3
clipLo(IC23v) <- -clipUp(IC23v)
stand(IC23v) <- as.matrix(stand3)
checkIC(IC23v)
Risks(IC23v) <- getFiRiskRegTS(risk = fiUnOvShoot(width = tau/sqrt(n)), 
                   ErrorDistr = Norm(), Regressor = K2, 
                   neighbor = TotalVarNeighborhood(radius = rad/2/sqrt(n)), 
                   clip = clipUp3, stand = stand3, sampleSize = n,
                   Algo = "A", cont = "left")
Infos(IC23v) <- matrix(c("manually", "O(n^(-1))-corrected solution"), ncol=2,
                    dimnames=list(character(0), c("method", "message")))
Risks(IC23v)
Risks(IC42v)
