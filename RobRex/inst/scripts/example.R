library(RobRex)
options("newDevice"=TRUE)

###############################################################################
## start of tests
###############################################################################

###############################################################################
## Example 1 (1-dim., discrete Regressor)
###############################################################################
K1 <- DiscreteDistribution(supp = 1:5)

# AL-estimator
system.time(IC.AL1 <- rgsOptIC.AL(r = 0.1, K = K1, check = TRUE), gcFirst = TRUE)
distrExOptions(ErelativeTolerance, 1e-10)
checkIC(IC.AL1)
distrExOptions(ErelativeTolerance, .Machine$double.eps^0.25)
Risks(IC.AL1)

# M-estimator
system.time(IC.M1 <- rgsOptIC.M(r = 0.1, K = K1, check = TRUE), gcFirst = TRUE)
checkIC(IC.M1)
Risks(IC.M1)

# MK-estimator
system.time(IC.MK1 <- rgsOptIC.MK(r = 0.1, K = K1, check = TRUE), gcFirst = TRUE)
checkIC(IC.MK1)
Risks(IC.MK1)

# ALc-estimator
system.time(IC.ALc1 <- rgsOptIC.ALc(r = 0.1, K = K1, check = TRUE), gcFirst = TRUE)
checkIC(IC.ALc1)
Risks(IC.ALc1)

# Mc-estimator
system.time(IC.Mc1 <- rgsOptIC.Mc(r = 0.1, K = K1, check = TRUE), gcFirst = TRUE)
checkIC(IC.Mc1)
Risks(IC.Mc1)

# ALs-estimator
system.time(IC.ALs1 <- rgsOptIC.ALs(r = 0.1, scale = 2, K = K1, check = TRUE), gcFirst = TRUE)
checkIC(IC.ALs1)
Risks(IC.ALs1)

# Ms-estimator
system.time(IC.Ms1 <- rgsOptIC.Ms(r = 0.1, K = K1, check = TRUE), gcFirst = TRUE)
checkIC(IC.Ms1)
Risks(IC.Ms1)

# BM-estimator
system.time(IC.BM1 <- rgsOptIC.BM(r = 0.1, K = K1), gcFirst = TRUE)
checkIC(IC.BM1)
Risks(IC.BM1)

# asymptotic MSEs
Risks(IC.AL1)$asMSE
Risks(IC.M1)$asMSE
Risks(IC.MK1)$asMSE
Risks(IC.ALc1)$asMSE
Risks(IC.Mc1)$asMSE
Risks(IC.ALs1)$asMSE
Risks(IC.Ms1)$asMSE
Risks(IC.BM1)$asMSE


###############################################################################
## Example 2 (1-dim., abs. cont. Regressor)
###############################################################################
K2 <- Unif(Min = 0, Max = 1)

# AL-estimator
IC.AL2 <- rgsOptIC.AL(r = 0.1, K = K2, check = TRUE)
#checkIC(IC.AL2, eval(CallL2Fam(IC.AL2)))
Risks(IC.AL2)

# M-estimator
IC.M2 <- rgsOptIC.M(r = 0.1, K = K2, check = TRUE)
#checkIC(IC.M2, eval(CallL2Fam(IC.M2)))
Risks(IC.M2)

# MK-estimator
IC.MK2 <- rgsOptIC.MK(r = 0.1, K = K2, check = TRUE)
#checkIC(IC.MK2, eval(CallL2Fam(IC.MK2)))
Risks(IC.MK2)

# ALc-estimator
# only implemented for discrete distributions

# Mc-estimator
# only implemented for 1-dim. discrete distributions

# ALs-estimator
IC.ALs2 <- rgsOptIC.ALs(r = 0.1, K = K2, check = TRUE)
#checkIC(IC.ALs2, eval(CallL2Fam(IC.ALs2)))
Risks(IC.ALs2)

# Ms-estimator
# only implemented for 1-dim. discrete distributions

# BM-estimator
# only defined for discrete distributions


# asymptotic MSEs
Risks(IC.AL2)$asMSE
Risks(IC.M2)$asMSE
Risks(IC.MK2)$asMSE
Risks(IC.ALs2)$asMSE


###############################################################################
## Example 3 (2-dim., discrete Regressor)
###############################################################################
K3 <- DiscreteMVDistribution(supp = matrix(c(0,1,1,1,0,2,2,0,2,1), ncol=2, byrow = TRUE))

# AL-estimator
IC.AL3 <- rgsOptIC.AL(r = 0.1, K = K3, check = TRUE)
checkIC(IC.AL3)
Risks(IC.AL3)

# M-estimator
IC.M3 <- rgsOptIC.M(r = 0.1, K = K3, check = TRUE)
checkIC(IC.M3)
Risks(IC.M3)

# MK-estimator
IC.MK3 <- rgsOptIC.MK(r = 0.1, K = K3, check = TRUE)
checkIC(IC.MK3)
Risks(IC.MK3)

# ALc-estimator
IC.ALc3 <- rgsOptIC.ALc(r = 0.1, K = K3, check = TRUE)
checkIC(IC.ALc3)
Risks(IC.ALc3)

# Mc-estimator
# only implemented for 1-dim. discrete distributions

# ALs-estimator
IC.ALs3 <- rgsOptIC.ALs(r = 0.1, K = K3, check = TRUE)
checkIC(IC.ALs3)
Risks(IC.ALs3)

# Ms-estimator
# only implemented for 1-dim. discrete distributions

# BM-estimator
# only defined for discrete distributions
# and only implemented for 1-dim. discrete distributions

# asymptotic MSEs
Risks(IC.AL3)$asMSE
Risks(IC.M3)$asMSE
Risks(IC.MK3)$asMSE
Risks(IC.ALc3)$asMSE
Risks(IC.ALs3)$asMSE


###############################################################################
## end of tests
###############################################################################

q("no")
