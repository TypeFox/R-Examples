###############################################################################
# Normal linear regression
###############################################################################
system.time(require(ROptRegTS))

tau <- qnorm(0.995)
n <- 2

# Regressor distributions
K1 <- DiscreteDistribution(c(0.5, 1.0, 1.5))
K2 <- Unif(Min = -1, Max = 2)

# generate a normal linear regression family
LM1 <- NormLinRegFamily(RegDistr = K1)
LM2 <- NormLinRegFamily(RegDistr = K2)

## classical optimal IC (radius = 0!)
RobLM1c <- InfRobRegTypeModel(center = LM1, 
                    neighbor = CondContNeighborhood(radiusCurve = function(x){rep(0, length(x))}))
RobLM1v <- InfRobRegTypeModel(center = LM1, 
                    neighbor = CondTotalVarNeighborhood(radiusCurve = function(x){rep(0, length(x))}))
RobLM2c <- InfRobRegTypeModel(center = LM2, 
                    neighbor = CondContNeighborhood(radiusCurve = function(x){rep(0, length(x))}))
RobLM2v <- InfRobRegTypeModel(center = LM2, 
                    neighbor = CondTotalVarNeighborhood(radiusCurve = function(x){rep(0, length(x))}))

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

###############################################################################
## asUnOvShoot solution
###############################################################################
## infinitesimal robust model
rad <- 0.4
RobLM1c <- InfRobRegTypeModel(center = LM1, 
                    neighbor = CondContNeighborhood(radiusCurve = function(x){rep(rad, length(x))}))
RobLM1v <- InfRobRegTypeModel(center = LM1, 
                    neighbor = CondTotalVarNeighborhood(radiusCurve = function(x){rep(rad/2, length(x))}))
RobLM2c <- InfRobRegTypeModel(center = LM2, 
                    neighbor = CondContNeighborhood(radiusCurve = function(x){rep(rad, length(x))}))
RobLM2v <- InfRobRegTypeModel(center = LM2, 
                    neighbor = CondTotalVarNeighborhood(radiusCurve = function(x){rep(rad/2, length(x))}))

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
RobLM3c <- FixRobRegTypeModel(center = LM1, 
                    neighbor = CondContNeighborhood(radiusCurve = function(x){rep(rad/sqrt(n), length(x))}))
RobLM3v <- FixRobRegTypeModel(center = LM1, 
                    neighbor = CondTotalVarNeighborhood(radiusCurve = function(x){rep(rad/2/sqrt(n), length(x))}))
RobLM4c <- FixRobRegTypeModel(center = LM2, 
                    neighbor = CondContNeighborhood(radiusCurve = function(x){rep(rad/sqrt(n), length(x))}))
RobLM4v <- FixRobRegTypeModel(center = LM2, 
                    neighbor = CondTotalVarNeighborhood(radiusCurve = function(x){rep(rad/2/sqrt(n), length(x))}))

## finite-sample optimal IC
system.time(IC32c <- optIC(model=RobLM3c, risk=fiUnOvShoot(width = tau/sqrt(n)), sampleSize = n), gcFirst = TRUE)
checkIC(IC32c) # takes some time
Risks(IC32c)

system.time(IC32v <- optIC(model=RobLM3v, risk=fiUnOvShoot(width = tau/sqrt(n)), sampleSize = n), gcFirst = TRUE)
checkIC(IC32v) # takes some time
Risks(IC32v)

distroptions("DefaultNrFFTGridPointsExponent", 8)
system.time(IC42c <- optIC(model=RobLM4c, risk=fiUnOvShoot(width = tau/sqrt(n)), sampleSize = n), gcFirst = TRUE)
checkIC(IC42c) # takes some time
Risks(IC42c)

system.time(IC42v <- optIC(model=RobLM4v, risk=fiUnOvShoot(width = tau/sqrt(n)), sampleSize = n), gcFirst = TRUE)
checkIC(IC42v) # takes some time
Risks(IC42v)

# function for computation of standardizing constant
Afct <- function(x, c0){return(x^2*pnorm(c0(x)))}

# O(n^(-0.5))-corrected solution 
# in case of contamination neighborhoods
IC13c <- IC12c
clipUp3 <- function(x){
    bf1 <- b1fun; eps <- radCurve
    b.as <- bf1(x)/(A*abs(x))
    return(pmax(0, b.as - eps(x)*(eps(x) + b.as*tau)/(sqrt(n)*2*tau*abs(x)*pnorm(-b.as))))
}
body(clipUp3) <- substitute({ bf1 <- b1fun; eps <- radCurve
                              b.as <- bf1(x)/(A*abs(x))
                              return(pmax(0, b.as - eps(x)*(eps(x) + b.as*tau)/(sqrt(n)*2*tau*abs(x)*pnorm(-b.as))))}, 
                            list(b1fun = clipUp(IC12c)@Map[[1]], 
                                 radCurve = function(x){rep(rad/sqrt(n), length(x))},
                                 A = as.vector(stand(IC12c)),
                                 tau = tau))
stand3 <- 1/(2*E(K1, Afct, c0 = clipUp3) - E(K1, function(x){x^2}))
clipUp13 <- function(x){ b3 <- b3fun; stand3*b3(x)*abs(x) }
body(clipUp13) <- substitute({ b3 <- b3fun; stand3*b3(x)*abs(x) }, list(b3fun = clipUp3, stand3 = stand3))
clipUp(IC13c) <- RealRandVariable(list(clipUp13), Domain=Reals())

clipLo13 <- function(x){ b3 <- b3fun; -stand3*b3(x)*abs(x) }
body(clipLo13) <- substitute({ b3 <- b3fun; -stand3*b3(x)*abs(x) }, list(b3fun = clipUp3, stand3 = stand3))
clipLo(IC13c) <- RealRandVariable(list(clipLo13), Domain=Reals())
stand(IC13c) <- as.matrix(stand3)
checkIC(IC13c)
distroptions("DefaultNrFFTGridPointsExponent", 12)
Risks(IC13c) <- getFiRiskRegTS(risk = fiUnOvShoot(width = tau/sqrt(n)), 
                   ErrorDistr = Norm(), Regressor = K1, 
                   neighbor = CondContNeighborhood(radius = 0, 
                                            radiusCurve = function(x){rep(rad/sqrt(n), length(x))}), 
                   clip = clipUp3, stand = stand3, sampleSize = n, cont = "left")
Risks(IC13c)
Infos(IC13c) <- matrix(c("manually", "O(n^(-1/2))-corrected solution"), ncol=2,
                    dimnames=list(character(0), c("method", "message")))
Risks(IC32c)

IC23c <- IC22c
clipUp3 <- function(x){
    bf1 <- b1fun; eps <- radCurve
    b.as <- bf1(x)/(A*abs(x))
    erg <- b.as - eps(x)*(eps(x) + b.as*tau)/(sqrt(n)*2*tau*abs(x)*pnorm(-b.as))
    return(pmax(0, erg))
}
body(clipUp3) <- substitute({ bf1 <- b1fun; eps <- radCurve
                              b.as <- bf1(x)/(A*abs(x))
                              erg <- b.as - eps(x)*(eps(x) + b.as*tau)/(sqrt(n)*2*tau*abs(x)*pnorm(-b.as))
                              return(pmax(0, erg))},
                            list(b1fun = clipUp(IC22c)@Map[[1]], 
                                 radCurve = function(x){rep(rad/sqrt(n), length(x))},
                                 A = as.vector(stand(IC22c)),
                                 tau = tau))
stand3 <- 1/(2*E(K2, Afct, c0 = clipUp3) - E(K2, function(x){x^2}))
clipUp23 <- function(x){ b3 <- b3fun; stand3*b3(x)*abs(x) }
body(clipUp23) <- substitute({ b3 <- b3fun; stand3*b3(x)*abs(x) }, list(b3fun = clipUp3, stand3 = stand3))
clipUp(IC23c) <- RealRandVariable(list(clipUp23), Domain=Reals())

clipLo23 <- function(x){ b3 <- b3fun; -stand3*b3(x)*abs(x) }
body(clipLo23) <- substitute({ b3 <- b3fun; -stand3*b3(x)*abs(x) }, list(b3fun = clipUp3, stand3 = stand3))
clipLo(IC23c) <- RealRandVariable(list(clipLo23), Domain=Reals())
stand(IC23c) <- as.matrix(stand3)
checkIC(IC23c)
distroptions("DefaultNrFFTGridPointsExponent", 8)
Risks(IC23c) <- getFiRiskRegTS(risk = fiUnOvShoot(width = tau/sqrt(n)), 
                   ErrorDistr = Norm(), Regressor = K2, 
                   neighbor = CondContNeighborhood(radius = 0, 
                                            radiusCurve = function(x){rep(rad/sqrt(n), length(x))}), 
                   clip = clipUp3, stand = stand3, sampleSize = n, cont = "left")
Risks(IC23c)
Infos(IC23c) <- matrix(c("manually", "O(n^(-1/2))-corrected solution"), ncol=2,
                    dimnames=list(character(0), c("method", "message")))
Risks(IC42c)

# O(n^(-1))-corrected solution 
# in case of total variation neighborhoods
IC13v <- IC12v
clipUp3 <- function(x){
    bf1 <- b1fun; delta <- radCurve
    b.as <- bf1(x)/(A*abs(x))
    return(pmax(0, b.as - tau*(2*b.as^2*delta(x) + tau*dnorm(b.as))/(n*6*pnorm(-b.as))))
}
body(clipUp3) <- substitute({ bf1 <- b1fun; delta <- radCurve
                              b.as <- bf1(x)/(A*abs(x))
                              return(pnorm(0, b.as - tau*(2*b.as^2*delta(x) + tau*dnorm(b.as))/(n*6*pnorm(-b.as))))}, 
                            list(b1fun = clipUp(IC12v)@Map[[1]], 
                                 radCurve = function(x){rep(rad/2/sqrt(n), length(x))},
                                 A = as.vector(stand(IC12v))))
stand3 <- 1/(2*E(K1, Afct, c0 = clipUp3) - E(K1, function(x){x^2}))
clipUp13 <- function(x){ b3 <- b3fun; stand3*b3(x)*abs(x) }
body(clipUp13) <- substitute({ b3 <- b3fun; stand3*b3(x)*abs(x) }, list(b3fun = clipUp3, stand3 = stand3))
clipUp(IC13v) <- RealRandVariable(list(clipUp13), Domain=Reals())

clipLo13 <- function(x){ b3 <- b3fun; -stand3*b3(x)*abs(x) }
body(clipLo13) <- substitute({ b3 <- b3fun; -stand3*b3(x)*abs(x) }, list(b3fun = clipUp3, stand3 = stand3))
clipLo(IC13v) <- RealRandVariable(list(clipLo13), Domain=Reals())
stand(IC13v) <- as.matrix(stand3)
checkIC(IC13v)

distroptions("DefaultNrFFTGridPointsExponent", 12)
Risks(IC13v) <- getFiRiskRegTS(risk = fiUnOvShoot(width = tau/sqrt(n)), 
                   ErrorDistr = Norm(), Regressor = K1, 
                   neighbor = CondTotalVarNeighborhood(radius = 0, 
                                            radiusCurve = function(x){rep(rad/2/sqrt(n), length(x))}), 
                   clip = clipUp3, stand = stand3, sampleSize = n, cont = "left")
Risks(IC13v)
Infos(IC13v) <- matrix(c("manually", "O(n^(-1))-corrected solution"), ncol=2,
                    dimnames=list(character(0), c("method", "message")))
Risks(IC32v)


IC23v <- IC22v
clipUp3 <- function(x){
    bf1 <- b1fun; delta <- radCurve
    b.as <- bf1(x)/(A*abs(x))
    return(pnorm(0, b.as - tau*(2*b.as^2*delta(x) + tau*dnorm(b.as))/(n*6*pnorm(-b.as))))
}
body(clipUp3) <- substitute({ bf1 <- b1fun; delta <- radCurve
                              b.as <- bf1(x)/(A*abs(x))
                              return(pnorm(0, b.as - tau*(2*b.as^2*delta(x) + tau*dnorm(b.as))/(n*6*pnorm(-b.as))))}, 
                            list(b1fun = clipUp(IC22v)@Map[[1]], 
                                 radCurve = function(x){rep(rad/2/sqrt(n), length(x))},
                                 A = as.vector(stand(IC22v))))
stand3 <- 1/(2*E(K2, Afct, c0 = clipUp3) - E(K2, function(x){x^2}))
clipUp23 <- function(x){ b3 <- b3fun; stand3*b3(x)*abs(x) }
body(clipUp23) <- substitute({ b3 <- b3fun; stand3*b3(x)*abs(x) }, list(b3fun = clipUp3, stand3 = stand3))
clipUp(IC23v) <- RealRandVariable(list(clipUp23), Domain=Reals())

clipLo23 <- function(x){ b3 <- b3fun; -stand3*b3(x)*abs(x) }
body(clipLo23) <- substitute({ b3 <- b3fun; -stand3*b3(x)*abs(x) }, list(b3fun = clipUp3, stand3 = stand3))
clipLo(IC23v) <- RealRandVariable(list(clipLo23), Domain=Reals())
stand(IC23v) <- as.matrix(stand3)
checkIC(IC23v)

distroptions("DefaultNrFFTGridPointsExponent", 8)
Risks(IC23v) <- getFiRiskRegTS(risk = fiUnOvShoot(width = tau/sqrt(n)), 
                   ErrorDistr = Norm(), Regressor = K2, 
                   neighbor = CondTotalVarNeighborhood(radius = 0, 
                                            radiusCurve = function(x){rep(rad/2/sqrt(n), length(x))}), 
                   clip = clipUp3, stand = stand3, sampleSize = n, cont = "left")
Risks(IC23v)
Infos(IC23v) <- matrix(c("manually", "O(n^(-1))-corrected solution"), ncol=2,
                    dimnames=list(character(0), c("method", "message")))
Risks(IC42v)
