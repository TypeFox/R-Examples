###############################################################################
# Normal linear regression
###############################################################################
require(ROptRegTS)

# Regressor distributions
K1 <- DiscreteDistribution(1:5)
K2 <- Unif(Min = 0, Max = 1)
K3 <- DiscreteMVDistribution(supp = matrix(c(0,1,1,1,0,2,2,0,2,1), ncol=2, byrow = TRUE))


# generate a normal linear regression family
(LM1 <- NormLinRegFamily(RegDistr = K1))
checkL2deriv(LM1)

(LM2 <- NormLinRegFamily(RegDistr = K2))
checkL2deriv(LM2)

(LM3 <- NormLinRegFamily(RegDistr = K3))
checkL2deriv(LM3)


###############################################################################
## classical optimal IC
###############################################################################
(IC01 <- optIC(model = LM1, risk = asCov()))
checkIC(IC01)
Risks(IC01)

(IC02 <- optIC(model = LM2, risk = asCov()))
checkIC(IC02)
Risks(IC02)

(IC03 <- optIC(model = LM3, risk = asCov()))
checkIC(IC03)
Risks(IC03)


## infinitesimal robust model
(RobLM1 <- InfRobRegTypeModel(center = LM1, neighbor = ContNeighborhood(radius = 0)))
(RobLM1c1 <- InfRobRegTypeModel(center = LM1, neighbor = Av1CondContNeighborhood(radius = 0)))
(RobLM1c2 <- InfRobRegTypeModel(center = LM1, neighbor = Av2CondContNeighborhood(radius = 0)))
(RobLM1v1 <- InfRobRegTypeModel(center = LM1, neighbor = Av1CondTotalVarNeighborhood(radius = 0)))

(RobLM2 <- InfRobRegTypeModel(center = LM2, neighbor = ContNeighborhood(radius = 0)))
(RobLM2c1 <- InfRobRegTypeModel(center = LM2, neighbor = Av1CondContNeighborhood(radius = 0)))
(RobLM2c2 <- InfRobRegTypeModel(center = LM2, neighbor = Av2CondContNeighborhood(radius = 0)))
(RobLM2v1 <- InfRobRegTypeModel(center = LM2, neighbor = Av1CondTotalVarNeighborhood(radius = 0)))

(RobLM3 <- InfRobRegTypeModel(center = LM3, neighbor = ContNeighborhood(radius = 0)))
(RobLM3c1 <- InfRobRegTypeModel(center = LM3, neighbor = Av1CondContNeighborhood(radius = 0)))
(RobLM3c2 <- InfRobRegTypeModel(center = LM3, neighbor = Av2CondContNeighborhood(radius = 0)))
(RobLM3v1 <- InfRobRegTypeModel(center = LM3, neighbor = Av1CondTotalVarNeighborhood(radius = 0)))


## classical optimal IC (radius = 0!)
system.time(IC11 <- optIC(model=RobLM1, risk=asMSE()), gcFirst = TRUE)
checkIC(IC11)
Risks(IC11)
system.time(IC11c1 <- optIC(model=RobLM1c1, risk=asMSE()), gcFirst = TRUE)
checkIC(IC11c1)
Risks(IC11c1)
system.time(IC11c2 <- optIC(model=RobLM1c2, risk=asMSE()), gcFirst = TRUE)
checkIC(IC11c2)
Risks(IC11c2)
system.time(IC11v1 <- optIC(model=RobLM1v1, risk=asMSE()), gcFirst = TRUE)
checkIC(IC11v1)
Risks(IC11v1)

system.time(IC12 <- optIC(model=RobLM2, risk=asMSE()), gcFirst = TRUE)
checkIC(IC12)
Risks(IC12)
system.time(IC12c1 <- optIC(model=RobLM2c1, risk=asMSE()), gcFirst = TRUE)
checkIC(IC12c1)
Risks(IC12c1)
system.time(IC12c2 <- optIC(model=RobLM2c2, risk=asMSE()), gcFirst = TRUE)
checkIC(IC12c2)
Risks(IC12c2)
system.time(IC12v1 <- optIC(model=RobLM2v1, risk=asMSE()), gcFirst = TRUE)
checkIC(IC12v1)
Risks(IC12v1)

system.time(IC13 <- optIC(model=RobLM3, risk=asMSE()), gcFirst = TRUE)
checkIC(IC13)
Risks(IC13)
system.time(IC13c1 <- optIC(model=RobLM3c1, risk=asMSE()), gcFirst = TRUE)
checkIC(IC13c1)
Risks(IC13c1)
system.time(IC13c2 <- optIC(model=RobLM3c2, risk=asMSE()), gcFirst = TRUE)
checkIC(IC13c2)
Risks(IC13c2)
system.time(IC13v1 <- optIC(model=RobLM3v1, risk=asMSE()), gcFirst = TRUE)
checkIC(IC13v1)
Risks(IC13v1)


###############################################################################
## MSE solution
###############################################################################
## infinitesimal robust model
(RobLM1 <- InfRobRegTypeModel(center = LM1, neighbor = ContNeighborhood(radius = 0.5)))
(RobLM1c1 <- InfRobRegTypeModel(center = LM1, neighbor = Av1CondContNeighborhood(radius = 0.5)))
(RobLM1c2 <- InfRobRegTypeModel(center = LM1, neighbor = Av2CondContNeighborhood(radius = 0.5)))
(RobLM1v1 <- InfRobRegTypeModel(center = LM1, neighbor = Av1CondTotalVarNeighborhood(radius = 0.25)))

(RobLM2 <- InfRobRegTypeModel(center = LM2, neighbor = ContNeighborhood(radius = 0.5)))
(RobLM2c1 <- InfRobRegTypeModel(center = LM2, neighbor = Av1CondContNeighborhood(radius = 0.5)))
(RobLM2c2 <- InfRobRegTypeModel(center = LM2, neighbor = Av2CondContNeighborhood(radius = 0.5)))
(RobLM2v1 <- InfRobRegTypeModel(center = LM2, neighbor = Av1CondTotalVarNeighborhood(radius = 0.25)))

(RobLM3 <- InfRobRegTypeModel(center = LM3, neighbor = ContNeighborhood(radius = 0.5)))
(RobLM3c1 <- InfRobRegTypeModel(center = LM3, neighbor = Av1CondContNeighborhood(radius = 0.5)))
(RobLM3c2 <- InfRobRegTypeModel(center = LM3, neighbor = Av2CondContNeighborhood(radius = 0.5)))
(RobLM3v1 <- InfRobRegTypeModel(center = LM3, neighbor = Av1CondTotalVarNeighborhood(radius = 0.25)))


## MSE solution
system.time(IC21 <- optIC(model=RobLM1, risk=asMSE()), gcFirst = TRUE)
checkIC(IC21)
Risks(IC21)
# without using symmetry
RobLM1@center@ErrorL2derivDistrSymm[[1]] <- NoSymmetry()
system.time(IC211 <- optIC(model=RobLM1, risk=asMSE()), gcFirst = TRUE)
checkIC(IC211)
Risks(IC211)
RobLM1@center@ErrorL2derivDistrSymm[[1]] <- SphericalSymmetry(SymmCenter = 0)

system.time(IC21c1 <- optIC(model=RobLM1c1, risk=asMSE()), gcFirst = TRUE)
checkIC(IC21c1)
Risks(IC21c1)
# without using symmetry
RobLM1c1@center@ErrorL2derivDistrSymm[[1]] <- NoSymmetry()
system.time(IC211c1 <- optIC(model=RobLM1c1, risk=asMSE()), gcFirst = TRUE)
checkIC(IC211c1)
Risks(IC211c1)
RobLM1c1@center@ErrorL2derivDistrSymm[[1]] <- SphericalSymmetry(SymmCenter = 0)

system.time(IC21c2 <- optIC(model=RobLM1c2, risk=asMSE()), gcFirst = TRUE)
checkIC(IC21c2)
Risks(IC21c2)
# without using symmetry
RobLM1c1@center@ErrorL2derivDistrSymm[[1]] <- NoSymmetry()
system.time(IC211c2 <- optIC(model=RobLM1c2, risk=asMSE()), gcFirst = TRUE)
checkIC(IC211c2)
Risks(IC211c2)
RobLM1c1@center@ErrorL2derivDistrSymm[[1]] <- SphericalSymmetry(SymmCenter = 0)

system.time(IC21v1 <- optIC(model=RobLM1v1, risk=asMSE()), gcFirst = TRUE)
checkIC(IC21v1)
Risks(IC21v1)
# without using symmetry
RobLM1v1@center@ErrorL2derivDistrSymm[[1]] <- NoSymmetry()
system.time(IC211v1 <- optIC(model=RobLM1v1, risk=asMSE()), gcFirst = TRUE)
checkIC(IC211v1)
Risks(IC211v1)
RobLM1v1@center@ErrorL2derivDistrSymm[[1]] <- SphericalSymmetry(SymmCenter = 0)

system.time(IC22 <- optIC(model=RobLM2, risk=asMSE()), gcFirst = TRUE)
checkIC(IC22)
Risks(IC22)
# without using symmetry
RobLM2@center@ErrorL2derivDistrSymm[[1]] <- NoSymmetry()
system.time(IC221 <- optIC(model=RobLM2, risk=asMSE()), gcFirst = TRUE)
checkIC(IC221)
Risks(IC221)
RobLM2@center@ErrorL2derivDistrSymm[[1]] <- SphericalSymmetry(SymmCenter = 0)

system.time(IC22c1 <- optIC(model=RobLM2c1, risk=asMSE()), gcFirst = TRUE)
checkIC(IC22c1)
Risks(IC22c1)
# without using symmetry
RobLM2c1@center@ErrorL2derivDistrSymm[[1]] <- NoSymmetry()
system.time(IC221c1 <- optIC(model=RobLM2c1, risk=asMSE()), gcFirst = TRUE)
checkIC(IC221c1)
Risks(IC221c1)
RobLM2c1@center@ErrorL2derivDistrSymm[[1]] <- SphericalSymmetry(SymmCenter = 0)

system.time(IC22c2 <- optIC(model=RobLM2c2, risk=asMSE()), gcFirst = TRUE)
checkIC(IC22c2)
Risks(IC22c2)
# without using symmetry
RobLM2c2@center@ErrorL2derivDistrSymm[[1]] <- NoSymmetry()
system.time(IC221c2 <- optIC(model=RobLM2c2, risk=asMSE()), gcFirst = TRUE)
checkIC(IC221c2)
Risks(IC221c2)
RobLM2c2@center@ErrorL2derivDistrSymm[[1]] <- SphericalSymmetry(SymmCenter = 0)

system.time(IC22v1 <- optIC(model=RobLM2v1, risk=asMSE()), gcFirst = TRUE)
checkIC(IC22v1)
Risks(IC22v1)
# without using symmetry
RobLM2v1@center@ErrorL2derivDistrSymm[[1]] <- NoSymmetry()
system.time(IC221v1 <- optIC(model=RobLM2v1, risk=asMSE()), gcFirst = TRUE)
checkIC(IC221v1)
Risks(IC221v1)
RobLM2v1@center@ErrorL2derivDistrSymm[[1]] <- SphericalSymmetry(SymmCenter = 0)

system.time(IC23 <- optIC(model=RobLM3, risk=asMSE()), gcFirst = TRUE)
checkIC(IC23)
Risks(IC23)
# without using symmetry
RobLM3@center@ErrorL2derivDistrSymm[[1]] <- NoSymmetry()
system.time(IC231 <- optIC(model=RobLM3, risk=asMSE()), gcFirst = TRUE)
checkIC(IC231)
Risks(IC231)
RobLM3@center@ErrorL2derivDistrSymm[[1]] <- SphericalSymmetry(SymmCenter = 0)

system.time(IC23c1 <- optIC(model=RobLM3c1, risk=asMSE()), gcFirst = TRUE)
checkIC(IC23c1)
Risks(IC23c1)
# without using symmetry
RobLM3c1@center@ErrorL2derivDistrSymm[[1]] <- NoSymmetry()
system.time(IC231c1 <- optIC(model=RobLM3c1, risk=asMSE()), gcFirst = TRUE)
checkIC(IC231c1)
Risks(IC231c1)
RobLM3c1@center@ErrorL2derivDistrSymm[[1]] <- SphericalSymmetry(SymmCenter = 0)

system.time(IC23c2 <- optIC(model=RobLM3c2, risk=asMSE()), gcFirst = TRUE)
checkIC(IC23c2)
Risks(IC23c2)
# without using symmetry
RobLM3c2@center@ErrorL2derivDistrSymm[[1]] <- NoSymmetry()
system.time(IC231c2 <- optIC(model=RobLM3c2, risk=asMSE()), gcFirst = TRUE)
checkIC(IC231c2)
Risks(IC231c2)
RobLM3c2@center@ErrorL2derivDistrSymm[[1]] <- SphericalSymmetry(SymmCenter = 0)

system.time(IC23v1 <- optIC(model=RobLM3v1, risk=asMSE()), gcFirst = TRUE)
checkIC(IC23v1)
Risks(IC23v1)
# without using symmetry
RobLM3v1@center@ErrorL2derivDistrSymm[[1]] <- NoSymmetry()
system.time(IC231v1 <- optIC(model=RobLM3v1, risk=asMSE()), gcFirst = TRUE)
checkIC(IC231v1)
Risks(IC231v1)
RobLM3v1@center@ErrorL2derivDistrSymm[[1]] <- SphericalSymmetry(SymmCenter = 0)


###############################################################################
## minimum bias solution
###############################################################################
system.time(IC31 <- optIC(model=RobLM1, risk=asBias(), tol = 1e-8), gcFirst = TRUE)
checkIC(IC31)
Risks(IC31)
# without using symmetry
RobLM1@center@ErrorL2derivDistrSymm[[1]] <- NoSymmetry()
system.time(IC311 <- optIC(model=RobLM1, risk=asBias(), tol = 1e-8), gcFirst = TRUE)
checkIC(IC311)
Risks(IC311)
RobLM1@center@ErrorL2derivDistrSymm[[1]] <- SphericalSymmetry(SymmCenter = 0)

system.time(IC31c1 <- optIC(model=RobLM1c1, risk=asBias(), tol = 1e-8), gcFirst = TRUE)
checkIC(IC31c1)
Risks(IC31c1)
# without using symmetry
RobLM1c1@center@ErrorL2derivDistrSymm[[1]] <- NoSymmetry()
system.time(IC311c1 <- optIC(model=RobLM1c1, risk=asBias(), tol = 1e-8), gcFirst = TRUE)
checkIC(IC311c1)
Risks(IC311c1)
RobLM1c1@center@ErrorL2derivDistrSymm[[1]] <- SphericalSymmetry(SymmCenter = 0)

system.time(IC31c2 <- optIC(model=RobLM1c2, risk=asBias(), tol = 1e-8), gcFirst = TRUE)
checkIC(IC31c2)
Risks(IC31c2)
# without using symmetry
RobLM1c2@center@ErrorL2derivDistrSymm[[1]] <- NoSymmetry()
system.time(IC311c2 <- optIC(model=RobLM1c2, risk=asBias(), tol = 1e-8), gcFirst = TRUE)
checkIC(IC311c2)
Risks(IC311c2)
RobLM1c2@center@ErrorL2derivDistrSymm[[1]] <- SphericalSymmetry(SymmCenter = 0)

system.time(IC31v1 <- optIC(model=RobLM1v1, risk=asBias(), tol = 1e-8), gcFirst = TRUE)
checkIC(IC31v1)
Risks(IC31v1)
# without using symmetry
RobLM1v1@center@ErrorL2derivDistrSymm[[1]] <- NoSymmetry()
system.time(IC311v1 <- optIC(model=RobLM1v1, risk=asBias(), tol = 1e-8), gcFirst = TRUE)
checkIC(IC311v1)
Risks(IC311v1)
RobLM1v1@center@ErrorL2derivDistrSymm[[1]] <- SphericalSymmetry(SymmCenter = 0)

system.time(IC32 <- optIC(model=RobLM2, risk=asBias(), tol = 1e-8), gcFirst = TRUE)
checkIC(IC32)
Risks(IC32)
# without using symmetry
RobLM2@center@ErrorL2derivDistrSymm[[1]] <- NoSymmetry()
system.time(IC321 <- optIC(model=RobLM2, risk=asBias(), tol = 1e-8), gcFirst = TRUE)
checkIC(IC321)
Risks(IC321)
RobLM2@center@ErrorL2derivDistrSymm[[1]] <- SphericalSymmetry(SymmCenter = 0)

system.time(IC32c1 <- optIC(model=RobLM2c1, risk=asBias(), tol = 1e-8), gcFirst = TRUE)
checkIC(IC32c1)
Risks(IC32c1)
# without using symmetry
RobLM2c1@center@ErrorL2derivDistrSymm[[1]] <- NoSymmetry()
system.time(IC321c1 <- optIC(model=RobLM2c1, risk=asBias(), tol = 1e-8), gcFirst = TRUE)
checkIC(IC321c1)
Risks(IC321c1)
RobLM2c1@center@ErrorL2derivDistrSymm[[1]] <- SphericalSymmetry(SymmCenter = 0)

system.time(IC32c2 <- optIC(model=RobLM2c2, risk=asBias(), tol = 1e-8), gcFirst = TRUE)
checkIC(IC32c2)
Risks(IC32c2)
# without using symmetry
RobLM2c2@center@ErrorL2derivDistrSymm[[1]] <- NoSymmetry()
system.time(IC321c2 <- optIC(model=RobLM2c2, risk=asBias(), tol = 1e-8), gcFirst = TRUE)
checkIC(IC321c2)
Risks(IC321c2)
RobLM2c2@center@ErrorL2derivDistrSymm[[1]] <- SphericalSymmetry(SymmCenter = 0)

system.time(IC32v1 <- optIC(model=RobLM2v1, risk=asBias(), tol = 1e-8), gcFirst = TRUE)
checkIC(IC32v1)
Risks(IC32v1)
# without using symmetry
RobLM2v1@center@ErrorL2derivDistrSymm[[1]] <- NoSymmetry()
system.time(IC321v1 <- optIC(model=RobLM2v1, risk=asBias(), tol = 1e-8), gcFirst = TRUE)
checkIC(IC321v1)
Risks(IC321v1)
RobLM2v1@center@ErrorL2derivDistrSymm[[1]] <- SphericalSymmetry(SymmCenter = 0)

system.time(IC33 <- optIC(model=RobLM3, risk=asBias(), tol = 1e-10), gcFirst = TRUE)
checkIC(IC33)
Risks(IC33)
# without using symmetry
RobLM3@center@ErrorL2derivDistrSymm[[1]] <- NoSymmetry()
system.time(IC331 <- optIC(model=RobLM3, risk=asBias(), tol = 1e-8), gcFirst = TRUE)
checkIC(IC331)
Risks(IC331)
RobLM3@center@ErrorL2derivDistrSymm[[1]] <- SphericalSymmetry(SymmCenter = 0)

system.time(IC33c1 <- optIC(model=RobLM3c1, risk=asBias(), tol = 1e-10), gcFirst = TRUE)
checkIC(IC33c1)
Risks(IC33c1)
# without using symmetry
RobLM3c1@center@ErrorL2derivDistrSymm[[1]] <- NoSymmetry()
system.time(IC331c1 <- optIC(model=RobLM3c1, risk=asBias(), tol = 1e-10), gcFirst = TRUE)
checkIC(IC331c1)
Risks(IC331c1)
RobLM3c1@center@ErrorL2derivDistrSymm[[1]] <- SphericalSymmetry(SymmCenter = 0)

system.time(IC33c2 <- optIC(model=RobLM3c2, risk=asBias(), tol = 1e-10), gcFirst = TRUE)
checkIC(IC33c2)
Risks(IC33c2)
# without using symmetry
RobLM3c2@center@ErrorL2derivDistrSymm[[1]] <- NoSymmetry()
system.time(IC331c2 <- optIC(model=RobLM3c2, risk=asBias(), tol = 1e-10), gcFirst = TRUE)
checkIC(IC331c2)
Risks(IC331c2)
RobLM3c2@center@ErrorL2derivDistrSymm[[1]] <- SphericalSymmetry(SymmCenter = 0)

system.time(IC33v1 <- optIC(model=RobLM3v1, risk=asBias(), tol = 1e-10), gcFirst = TRUE)
checkIC(IC33v1)
Risks(IC33v1)
# without using symmetry
RobLM3v1@center@ErrorL2derivDistrSymm[[1]] <- NoSymmetry()
system.time(IC331v1 <- optIC(model=RobLM3v1, risk=asBias(), tol = 1e-10), gcFirst = TRUE)
checkIC(IC331v1)
Risks(IC331v1)
RobLM3v1@center@ErrorL2derivDistrSymm[[1]] <- SphericalSymmetry(SymmCenter = 0)


###############################################################################
## radius minimax IC
###############################################################################
system.time(IC41 <- radiusMinimaxIC(L2Fam = LM1, neighbor = ContNeighborhood(), 
                                    risk = asMSE(), loRad = 0, upRad = 2.0), gcFirst = TRUE)
checkIC(IC41)
Risks(IC41)
system.time(IC41c1 <- radiusMinimaxIC(L2Fam = LM1, neighbor = Av1CondContNeighborhood(), 
                                    risk = asMSE(), loRad = 0, upRad = Inf), gcFirst = TRUE)
checkIC(IC41c1)
Risks(IC41c1)
system.time(IC41c2 <- radiusMinimaxIC(L2Fam = LM1, neighbor = Av2CondContNeighborhood(), 
                                    risk = asMSE(), loRad = 0, upRad = Inf), gcFirst = TRUE)
checkIC(IC41c2)
Risks(IC41c2)
system.time(IC41v1 <- radiusMinimaxIC(L2Fam = LM1, neighbor = Av1CondTotalVarNeighborhood(), 
                                    risk = asMSE(), loRad = 0, upRad = 2), gcFirst = TRUE)
checkIC(IC41v1)
Risks(IC41v1)

system.time(IC42 <- radiusMinimaxIC(L2Fam = LM2, neighbor = ContNeighborhood(), 
                                    risk = asMSE(), loRad = 0, upRad = Inf), gcFirst = TRUE)
checkIC(IC42)
Risks(IC42)
system.time(IC42c1 <- radiusMinimaxIC(L2Fam = LM2, neighbor = Av1CondContNeighborhood(), 
                                    risk = asMSE(), loRad = 0, upRad = 1.0), gcFirst = TRUE)
checkIC(IC42c1)
Risks(IC42c1)
system.time(IC42c2 <- radiusMinimaxIC(L2Fam = LM2, neighbor = Av2CondContNeighborhood(), 
                                    risk = asMSE(), loRad = 0, upRad = 1.0), gcFirst = TRUE)
checkIC(IC42c2)
Risks(IC42c2)
system.time(IC42v1 <- radiusMinimaxIC(L2Fam = LM2, neighbor = Av1CondTotalVarNeighborhood(), 
                                    risk = asMSE(), loRad = 0, upRad = 1.0), gcFirst = TRUE)
checkIC(IC42v1)
Risks(IC42v1)

system.time(IC43 <- radiusMinimaxIC(L2Fam = LM3, neighbor = ContNeighborhood(), 
                                    risk = asMSE(), loRad = 0, upRad = Inf), gcFirst = TRUE)
checkIC(IC43)
Risks(IC43)
system.time(IC43c1 <- radiusMinimaxIC(L2Fam = LM3, neighbor = Av1CondContNeighborhood(), 
                                    risk = asMSE(), loRad = 0, upRad = 1.0), gcFirst = TRUE)
checkIC(IC43c1)
Risks(IC43c1)
system.time(IC43c2 <- radiusMinimaxIC(L2Fam = LM3, neighbor = Av2CondContNeighborhood(), 
                                    risk = asMSE(), loRad = 0, upRad = Inf), gcFirst = TRUE)
checkIC(IC43c2)
Risks(IC43c2)
system.time(IC43v1 <- radiusMinimaxIC(L2Fam = LM3, neighbor = Av1CondTotalVarNeighborhood(), 
                                    risk = asMSE(), loRad = 0, upRad = 1.0), gcFirst = TRUE)
checkIC(IC43v1)
Risks(IC43v1)


###############################################################################
## least favorable radius
###############################################################################
system.time(r.rho1 <- leastFavorableRadius(L2Fam = LM1, neighbor = ContNeighborhood(), 
                                           risk = asMSE(), rho = 0.5), gcFirst = TRUE)
r.rho1
system.time(r.rho1c1 <- leastFavorableRadius(L2Fam = LM1, neighbor = Av1CondContNeighborhood(), 
                                            risk = asMSE(), rho = 1/3), gcFirst = TRUE)
r.rho1c1
# compare 1-dim. normal location
system.time(r.rho1c2 <- leastFavorableRadius(L2Fam = LM1, neighbor = Av2CondContNeighborhood(), 
                                            risk = asMSE(), rho = 1/3), gcFirst = TRUE)
r.rho1c2
system.time(r.rho1v1 <- leastFavorableRadius(L2Fam = LM1, neighbor = Av1CondTotalVarNeighborhood(), 
                                            risk = asMSE(), rho = 1/3), gcFirst = TRUE)
r.rho1v1

system.time(r.rho2 <- leastFavorableRadius(L2Fam = LM2, neighbor = ContNeighborhood(), 
                                           risk = asMSE(), rho = 0.25), gcFirst = TRUE)
r.rho2
system.time(r.rho2c1 <- leastFavorableRadius(L2Fam = LM2, neighbor = Av1CondContNeighborhood(), 
                                            risk = asMSE(), rho = 0.75), gcFirst = TRUE)
r.rho2c1
# compare 1-dim. normal location
system.time(r.rho2c2 <- leastFavorableRadius(L2Fam = LM2, neighbor = Av2CondContNeighborhood(), 
                                            risk = asMSE(), rho = 0.5), gcFirst = TRUE)
r.rho2c2
system.time(r.rho2v1 <- leastFavorableRadius(L2Fam = LM2, neighbor = Av1CondTotalVarNeighborhood(), 
                                            risk = asMSE(), rho = 0.75), gcFirst = TRUE)
r.rho2v1

system.time(r.rho3 <- leastFavorableRadius(L2Fam = LM3, neighbor = ContNeighborhood(), 
                                           risk = asMSE(), rho = 0.25), gcFirst = TRUE)
r.rho3
system.time(r.rho3c1 <- leastFavorableRadius(L2Fam = LM3, neighbor = Av1CondContNeighborhood(), 
                                            risk = asMSE(), rho = 0.5), gcFirst = TRUE)
r.rho3c1
# compare 1-dim. normal location
system.time(r.rho3c2 <- leastFavorableRadius(L2Fam = LM3, neighbor = Av2CondContNeighborhood(), 
                                            risk = asMSE(), rho = 1/3), gcFirst = TRUE)
r.rho3c2
system.time(r.rho3v1 <- leastFavorableRadius(L2Fam = LM3, neighbor = Av1CondTotalVarNeighborhood(), 
                                            risk = asMSE(), rho = 0.5), gcFirst = TRUE)
r.rho3v1
