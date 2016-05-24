###############################################################################
# Normal linear regression with unknown scale
###############################################################################

require(ROptRegTS)

###############################################################################
## Example 1
###############################################################################
a0 <- 0.25
K1 <- DiscreteDistribution(supp = c(-1, 1+a0), prob = c((1+a0)/(2+a0), 1/(2+a0)))

(LM11 <- NormLinRegFamily(RegDistr = K1))
checkL2deriv(LM11)

(LM12 <- NormLinRegInterceptFamily(RegDistr = K1, trafo = matrix(c(1,0), nrow = 1), nuisance = TRUE))
checkL2deriv(LM12)

## infinitesimal robust model
(RobLM11 <- InfRobRegTypeModel(center = LM11, neighbor = ContNeighborhood(radius = 0.5)))
(RobLM12 <- InfRobRegTypeModel(center = LM12, neighbor = ContNeighborhood(radius = 0.5)))

## MSE solution
system.time(IC11 <- optIC(model=RobLM11, risk=asMSE()), gcFirst = TRUE)
checkIC(IC11)
Risks(IC11)

system.time(IC12 <- optIC(model=RobLM12, risk=asMSE(), tol = 1e-4), gcFirst = TRUE)
checkIC(IC12)
Risks(IC12)
(M1 <- -stand(IC12)[,2]/stand(IC12)[,1])


###############################################################################
## Example 2 
###############################################################################
a0 <- 0.25
K2 <- DiscreteDistribution(supp = c(-1, 0, 1+a0), prob = c(1/(2+a0), a0/(1+a0), 1/((1+a0)*(2+a0))))

(LM21 <- NormLinRegFamily(RegDistr = K2))
checkL2deriv(LM21)

(LM22 <- NormLinRegInterceptFamily(RegDistr = K2, trafo = matrix(c(1,0), nrow = 1), nuisance = TRUE))
checkL2deriv(LM22)

## infinitesimal robust model
(RobLM21 <- InfRobRegTypeModel(center = LM21, neighbor = ContNeighborhood(radius = 0.5)))
(RobLM22 <- InfRobRegTypeModel(center = LM22, neighbor = ContNeighborhood(radius = 0.5)))

## MSE solution
system.time(IC21 <- optIC(model=RobLM21, risk=asMSE()), gcFirst = TRUE)
checkIC(IC21)
Risks(IC21)

system.time(IC22 <- optIC(model=RobLM22, risk=asMSE()), gcFirst = TRUE)
checkIC(IC22)
Risks(IC22)
(M2 <- -stand(IC22)[,2]/stand(IC22)[,1])


###############################################################################
## Example 3 
###############################################################################
p <- 0.55
K3 <- ConvexContamination(e1 = Unif(Min = -1, Max = 0), 
                          e2 = Unif(Min = 0, Max = p/(1-p)), size = 1-p)

(LM31 <- NormLinRegFamily(RegDistr = K3))
checkL2deriv(LM31)

distrExOptions(ErelativeTolerance, .Machine$double.eps^0.5)
(LM32 <- NormLinRegInterceptFamily(RegDistr = K3, trafo = matrix(c(1,0), nrow = 1), nuisance = TRUE))
distrExOptions(ErelativeTolerance, .Machine$double.eps^0.25)
checkL2deriv(LM32)

## infinitesimal robust model
(RobLM31 <- InfRobRegTypeModel(center = LM31, neighbor = ContNeighborhood(radius = 0.5)))
(RobLM32 <- InfRobRegTypeModel(center = LM32, neighbor = ContNeighborhood(radius = 0.5)))

## MSE solution
system.time(IC31 <- optIC(model=RobLM31, risk=asMSE()), gcFirst = TRUE)
checkIC(IC31)
Risks(IC31)

# takes quite some time!
system.time(IC32 <- optIC(model=RobLM32, risk=asMSE()), gcFirst = TRUE)
checkIC(IC32)
Risks(IC32)
(M3 <- -stand(IC32)[,2]/stand(IC32)[,1])
