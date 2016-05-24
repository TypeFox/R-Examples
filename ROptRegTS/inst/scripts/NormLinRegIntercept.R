###############################################################################
# Normal linear regression with unknown intercept
###############################################################################

require(ROptRegTS)

###############################################################################
## Example 1 (1-dim., discrete Regressor)
###############################################################################
K1 <- DiscreteDistribution(supp = c(-2, -1, 1, 2))

(LM1 <- NormLinRegInterceptFamily(RegDistr = K1))
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

system.time(IC21c <- optIC(model=RobLM1c, risk=asMSE()), gcFirst = TRUE)
checkIC(IC21c)
Risks(IC21c)

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
K2 <- Unif(Min = -1, Max = 1)

(LM2 <- NormLinRegInterceptFamily(RegDistr = K2))
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

## MSE solution
system.time(IC22c <- optIC(model=RobLM2c, risk=asMSE()), gcFirst = TRUE)
checkIC(IC22c)
Risks(IC22c)

## minimum bias solution
system.time(IC32 <- optIC(model=RobLM2, risk=asBias(), tol = 1e-10), gcFirst = TRUE)
checkIC(IC32)
Risks(IC32)
