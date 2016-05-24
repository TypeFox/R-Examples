###############################################################################
# Normal linear regression with unknown scale
###############################################################################

require(ROptRegTS)
require(RobRex)

###############################################################################
## Example 1 (1-dim., discrete Regressor)
###############################################################################
K1 <- DiscreteDistribution(supp = 1:5)

(LM1 <- NormLinRegScaleFamily(RegDistr = K1))
checkL2deriv(LM1)

(IC01 <- optIC(model = LM1, risk = asCov()))
checkIC(IC01)
Risks(IC01)

## infinitesimal robust model
(RobLM1 <- InfRobRegTypeModel(center = LM1, neighbor = ContNeighborhood(radius = 0)))
(RobLM1c <- InfRobRegTypeModel(center = LM1, neighbor = Av1CondContNeighborhood(radius = 0)))

## MSE solution
system.time(IC11 <- optIC(model=RobLM1, risk=asMSE()), gcFirst = TRUE)
checkIC(IC11)
Risks(IC11)
system.time(IC11c <- optIC(model=RobLM1c, risk=asMSE()), gcFirst = TRUE)
checkIC(IC11c)
Risks(IC11c)

## infinitesimal robust model
(RobLM1 <- InfRobRegTypeModel(center = LM1, neighbor = ContNeighborhood(radius = 0.5)))
(RobLM1c <- InfRobRegTypeModel(center = LM1, neighbor = Av1CondContNeighborhood(radius = 0.5)))

## MSE solution
system.time(IC21 <- optIC(model=RobLM1, risk=asMSE()), gcFirst = TRUE)
checkIC(IC21)
Risks(IC21)

# AL-estimator from package RobRex
system.time(IC.AL1 <- rgsOptIC.AL(r = 0.5, K = K1, check = TRUE), gcFirst = TRUE)
checkIC(IC.AL1)
Risks(IC.AL1)

## MSE solution
system.time(IC21c <- optIC(model=RobLM1c, risk=asMSE()), gcFirst = TRUE)
checkIC(IC21c)
Risks(IC21c)

# AL-estimator from package RobRex
system.time(IC.AL1c <- rgsOptIC.ALc(r = 0.5, K = K1, check = TRUE), gcFirst = TRUE)
checkIC(IC.AL1c)
Risks(IC.AL1c)

## minimum bias solution
system.time(IC31 <- optIC(model=RobLM1, risk=asBias(), tol = 1e-10), gcFirst = TRUE)
checkIC(IC31)
Risks(IC31)

## radius minimax IC (takes quite some time!)
system.time(IC41 <- radiusMinimaxIC(L2Fam = LM1, neighbor = ContNeighborhood(), 
                                    risk = asMSE(), loRad = 0, upRad = 2.0), gcFirst = TRUE)
checkIC(IC41)
Risks(IC41)
system.time(IC41c <- radiusMinimaxIC(L2Fam = LM1, neighbor = Av1CondContNeighborhood(), 
                                    risk = asMSE(), loRad = 0, upRad = 2.0), gcFirst = TRUE)
checkIC(IC41c)
Risks(IC41c)

## least favorable radius (takes a very long time!)
#system.time(r.rho <- leastFavorableRadius(L2Fam = LM1, neighbor = ContNeighborhood(), 
#                                          risk = asMSE(), rho = 0.5), gcFirst = TRUE)
#r.rho
#system.time(r.rhoc <- leastFavorableRadius(L2Fam = LM1, neighbor = Av1CondContNeighborhood(), 
#                                            risk = asMSE(), rho = 1/3), gcFirst = TRUE)
#r.rhoc



###############################################################################
## Example 2 (1-dim., abs. cont. Regressor)
###############################################################################
K2 <- Unif(Min = 0, Max = 1)

(LM2 <- NormLinRegScaleFamily(RegDistr = K2))
checkL2deriv(LM2)

(IC02 <- optIC(model = LM2, risk = asCov()))
checkIC(IC02)
Risks(IC02)

## infinitesimal robust model
(RobLM2 <- InfRobRegTypeModel(center = LM2, neighbor = ContNeighborhood(radius = 0)))
(RobLM2c <- InfRobRegTypeModel(center = LM2, neighbor = Av1CondContNeighborhood(radius = 0)))

## MSE solution
system.time(IC12 <- optIC(model=RobLM2, risk=asMSE()), gcFirst = TRUE)
checkIC(IC12)
Risks(IC12)
system.time(IC12c <- optIC(model=RobLM2c, risk=asMSE()), gcFirst = TRUE)
checkIC(IC12c)
Risks(IC12c)

## infinitesimal robust model
(RobLM2 <- InfRobRegTypeModel(center = LM2, neighbor = ContNeighborhood(radius = 0.5)))
(RobLM2c <- InfRobRegTypeModel(center = LM2, neighbor = Av1CondContNeighborhood(radius = 0.5)))

## MSE solution
system.time(IC22 <- optIC(model=RobLM2, risk=asMSE()), gcFirst = TRUE)
checkIC(IC22)
Risks(IC22)

# AL-estimator from package RobRex
system.time(IC.AL2 <- rgsOptIC.AL(r = 0.5, K = K2, check = TRUE), gcFirst = TRUE)
checkIC(IC.AL2)
Risks(IC.AL2)

## MSE solution
system.time(IC22c <- optIC(model=RobLM2c, risk=asMSE()), gcFirst = TRUE)
checkIC(IC22c)
Risks(IC22c)


## MSE solution
system.time(IC32 <- optIC(model=RobLM2, risk=asBias(), tol = 1e-10), gcFirst = TRUE)
checkIC(IC32)
Risks(IC32)


###############################################################################
## Example 3 (2-dim., discrete Regressor)
###############################################################################
K3 <- DiscreteMVDistribution(supp = matrix(c(0,1,1,1,0,2,2,0,2,1), ncol=2, byrow = TRUE))

(LM3 <- NormLinRegScaleFamily(RegDistr = K3))
checkL2deriv(LM3)

(IC03 <- optIC(model = LM3, risk = asCov()))
checkIC(IC03)
Risks(IC03)

## infinitesimal robust model
(RobLM3 <- InfRobRegTypeModel(center = LM3, neighbor = ContNeighborhood(radius = 0)))
(RobLM3c <- InfRobRegTypeModel(center = LM3, neighbor = Av1CondContNeighborhood(radius = 0)))

## MSE solution
system.time(IC13 <- optIC(model=RobLM3, risk=asMSE()), gcFirst = TRUE)
checkIC(IC13)
Risks(IC13)
system.time(IC13c <- optIC(model=RobLM3c, risk=asMSE()), gcFirst = TRUE)
checkIC(IC13c)
Risks(IC13c)

## infinitesimal robust model
(RobLM3 <- InfRobRegTypeModel(center = LM3, neighbor = ContNeighborhood(radius = 0.5)))
(RobLM3c <- InfRobRegTypeModel(center = LM3, neighbor = Av1CondContNeighborhood(radius = 0.5)))

## MSE solution
system.time(IC23 <- optIC(model=RobLM3, risk=asMSE()), gcFirst = TRUE)
checkIC(IC23)
Risks(IC23)

# AL-estimator from package RobRex
system.time(IC.AL3 <- rgsOptIC.AL(r = 0.5, K = K3, check = TRUE), gcFirst = TRUE)
checkIC(IC.AL3)
Risks(IC.AL3)

## MSE solution
system.time(IC23c <- optIC(model=RobLM3c, risk=asMSE(), tol=5e-6), gcFirst = TRUE)
checkIC(IC23c)
Risks(IC23c)

# AL-estimator from package RobRex
system.time(IC.AL3c <- rgsOptIC.ALc(r = 0.5, K = K3, check = TRUE), gcFirst = TRUE)
checkIC(IC.AL3c)
Risks(IC.AL3c)

## minimum bias solution
system.time(IC33 <- optIC(model=RobLM3, risk=asBias(), tol = 1e-10), gcFirst = TRUE)
checkIC(IC33)
Risks(IC33)
