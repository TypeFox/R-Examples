# ==============================================================================
# STOICHCALC - R-Routines for Solving Stoichiometric Equations
# ==============================================================================
#
# Realistic example                                                   22.02.2010
# =================
#
# Literature: Peter Reichert and Nele Schuwirth
#             A generic framework for deriving process stoichiometry 
#             in environmental models
#             submitted to Environmental Modelling and Software
#
# ==============================================================================


# Load library functions:
# =======================

source("stoichcalc.r")


# Define composition:
# ===================

# define composition parameters:
# ------------------------------

param <- list(alpha.O.ALG     = 0.50,         # gO/gALG
              alpha.H.ALG     = 0.07,         # gH/gALG
              alpha.N.ALG     = 0.06,         # gN/gALG
              alpha.P.ALG     = 0.005,        # gP/gALG
              alpha.O.ZOO     = 0.50,         # gO/gZOO
              alpha.H.ZOO     = 0.07,         # gH/gZOO
              alpha.N.ZOO     = 0.06,         # gN/gZOO
              alpha.P.ZOO     = 0.01,         # gP/gZOO
              alpha.O.POM     = 0.40,         # gO/gPOM
              alpha.H.POM     = 0.07,         # gH/gPOM
              alpha.N.POM     = 0.04,         # gN/gPOM
              alpha.P.POM     = 0.007,        # gP/gPOM
              alpha.O.DOM     = 0.40,         # gO/gDOM
              alpha.H.DOM     = 0.07,         # gH/gDOM
              alpha.N.DOM     = 0.04,         # gN/gDOM
              alpha.P.DOM     = 0.007,        # gP/gDOM
              Y.ZOO           = 0.2,          # gZOO/gALG
              f.DOM           = 0.1,          # gDOM/gALG
              f.POM           = 0.2)          # gPOM/gALG
            
# calculate carbon fractions to guarantee that the fractions sum to unity:

param$alpha.C.ALG <- 1 - (param$alpha.O.ALG + param$alpha.H.ALG + 
                          param$alpha.N.ALG + param$alpha.P.ALG)
param$alpha.C.ZOO <- 1 - (param$alpha.O.ZOO + param$alpha.H.ZOO + 
                          param$alpha.N.ZOO + param$alpha.P.ZOO)
param$alpha.C.POM <- 1 - (param$alpha.O.POM + param$alpha.H.POM + 
                          param$alpha.N.POM + param$alpha.P.POM)
param$alpha.C.DOM <- 1 - (param$alpha.O.DOM + param$alpha.H.DOM + 
                          param$alpha.N.DOM + param$alpha.P.DOM)


# define composition vectors of substances an organisms:
# ------------------------------------------------------

NH4    <- c(H      = 4*1/14,              # gH/gNH4-N
            N      = 1,                   # gN/gNH4-N
            charge = 1/14)                # chargeunits/gNH4-N
NO3    <- c(O      = 3*16/14,             # gO/gNO3-N
            N      = 1,                   # gN/gNO3-N
            charge = -1/14)               # chargeunits/gNO3-N
HPO4   <- c(O      = 4*16/31,             # gO/gHPO4-P
            H      = 1*1/31,              # gH/gHPO4-P
            P      = 1,                   # gP/gHPO4-P
            charge = -2/31)               # chargeunits/gHPO4-P
HCO3   <- c(C      = 1,                   # gC/gHCO3-C
            O      = 3*16/12,             # gO/gHCO3-C
            H      = 1*1/12,              # gH/gHCO3-C
            charge = -1/12)               # chargeunits/gHCO3-C
O2     <- c(O      = 1)                   # gO/gO2-O
H      <- c(H      = 1,                   # gH/molH
            charge = 1)                   # chargeunits/molH
H2O    <- c(O      = 1*12,                # gO/molH2O
            H      = 2*1)                 # gH/molH2O
ALG    <- c(C      = param$alpha.C.ALG,     # gC/gALG
            O      = param$alpha.O.ALG,     # gO/gALG
            H      = param$alpha.H.ALG,     # gH/gALG
            N      = param$alpha.N.ALG,     # gN/gALG
            P      = param$alpha.P.ALG)     # gP/gALG
ZOO    <- c(C      = param$alpha.C.ZOO,     # gC/gZOO
            O      = param$alpha.O.ZOO,     # gO/gZOO
            H      = param$alpha.H.ZOO,     # gH/gZOO
            N      = param$alpha.N.ZOO,     # gN/gZOO
            P      = param$alpha.P.ZOO)     # gP/gZOO
POM    <- c(C      = param$alpha.C.POM,     # gC/gPOM
            O      = param$alpha.O.POM,     # gO/gPOM
            H      = param$alpha.H.POM,     # gH/gPOM
            N      = param$alpha.N.POM,     # gN/gPOM
            P      = param$alpha.P.POM)     # gP/gPOM
DOM    <- c(C      = param$alpha.C.DOM,     # gC/gDOM
            O      = param$alpha.O.DOM,     # gO/gDOM
            H      = param$alpha.H.DOM,     # gH/gDOM
            N      = param$alpha.N.DOM,     # gN/gDOM
            P      = param$alpha.P.DOM)     # gP/gDOM


# define list of composition vectors:
# -----------------------------------

subst.comp <- list(NH4  = NH4,
                   NO3  = NO3,
                   HPO4 = HPO4,
                   HCO3 = HCO3,
                   O2   = O2,
                   H    = H,
                   H2O  = H2O,
                   ALG  = ALG,
                   ZOO  = ZOO,
                   POM  = POM,
                   DOM  = DOM)


# compile and print composition matrix:
# -------------------------------------
 
alpha <- calc.comp.matrix(subst.comp)

print(alpha)

nu.basis <- calc.stoich.basis(alpha)

print(nu.basis)

#test
nu.basis %*% t(alpha)


# Derivation of Process Stoichiometry:
# ====================================

# test of error messages for growth of zooplankton: 
# -------------------------------------------------

# substances/organisms relevant for growth of zooplankton:

subst.gro.ZOO <- c("NH4","HPO4","HCO3","O2","H","H2O","ALG","ZOO","POM","DOM")

# omission of constraints:

basis.gro.ZOO <- calc.stoich.basis(alpha       = alpha,
                                   subst       = subst.gro.ZOO)
                                   
nu.gro.ZOO    <- calc.stoich.coef (alpha       = alpha,
                                   name        = "gro.ZOO",
                                   subst       = subst.gro.ZOO,
                                   subst.norm  = "ZOO",
                                   nu.norm     = 1)
                                
# insufficient number of constraints:

const.gro.ZOO <- list(c("ZOO" = 1,"ALG" = param$Y.ZOO)) 

basis.gro.ZOO <- calc.stoich.basis(alpha       = alpha,
                                   subst       = subst.gro.ZOO,
                                   constraints = const.gro.ZOO)
                           
nu.gro.ZOO    <- calc.stoich.coef (alpha       = alpha,
                                   name        = "gro.ZOO",
                                   subst       = subst.gro.ZOO,
                                   subst.norm  = "ZOO",
                                   nu.norm     = 1,
                                   constraints = const.gro.ZOO)

# wrong substance name for normalization: 

const.gro.ZOO <- list(c("ZOO" = 1,"ALG" = param$Y.ZOO),
                      c("POM" = 1,"ALG" = param$f.POM),
                      c("DOM" = 1,"ALG" = param$f.DOM)) 

basis.gro.ZOO <- calc.stoich.basis(alpha       = alpha,
                                   subst       = subst.gro.ZOO,
                                   constraints = const.gro.ZOO)
                           
nu.gro.ZOO    <- calc.stoich.coef (alpha       = alpha,
                                   name        = "gro.ZOO",
                                   subst       = subst.gro.ZOO,
                                   subst.norm  = "FOO",
                                   nu.norm     = 1,
                                   constraints = constraints)

# good process definition 

const.gro.ZOO <- list(c("ZOO" = 1,"ALG" = param$Y.ZOO),
                      c("POM" = 1,"ALG" = param$f.POM),
                      c("DOM" = 1,"ALG" = param$f.DOM)) 

basis.gro.ZOO <- calc.stoich.basis(alpha       = alpha,
                                   subst       = subst.gro.ZOO,
                                   constraints = const.gro.ZOO)
                           
nu.gro.ZOO    <- calc.stoich.coef (alpha       = alpha,
                                   name        = "gro.ZOO",
                                   subst       = subst.gro.ZOO,
                                   subst.norm  = "ZOO",
                                   nu.norm     = 1,
                                   constraints = const.gro.ZOO)
                                   
print(nu.gro.ZOO)

# demonstration that stoichiometric coefficient of PO4 becomes negative
# if the P content of phytoplankton is too low:

alpha.test <- alpha
alpha.test["P","ALG"] <- 0.004

const.gro.ZOO <- list(c("ZOO" = 1,"ALG" = param$Y.ZOO),
                      c("POM" = 1,"ALG" = param$f.POM),
                      c("DOM" = 1,"ALG" = param$f.DOM)) 

nu.gro.ZOO    <- calc.stoich.coef (alpha       = alpha.test,
                                   name        = "gro.ZOO",
                                   subst       = subst.gro.ZOO,
                                   subst.norm  = "ZOO",
                                   nu.norm     = 1,
                                   constraints = const.gro.ZOO)

print(nu.gro.ZOO)

# this can be dorrected by making the yield smaller:

const.gro.ZOO <- list(c("ZOO" = 1,"ALG" = param$Y.ZOO/2),
                      c("POM" = 1,"ALG" = param$f.POM),
                      c("DOM" = 1,"ALG" = param$f.DOM)) 

nu.gro.ZOO    <- calc.stoich.coef (alpha       = alpha.test,
                                   name        = "gro.ZOO",
                                   subst       = subst.gro.ZOO,
                                   subst.norm  = "ZOO",
                                   nu.norm     = 1,
                                   constraints = const.gro.ZOO)

print(nu.gro.ZOO)

# calculating stoichiometric matrix: 
# ----------------------------------

# growth of zooplankton:

subst.gro.ZOO <- c("NH4","HPO4","HCO3","O2","H","H2O","ALG","ZOO","POM","DOM")

const.gro.ZOO <- list(c("ZOO" = 1,"ALG" = param$Y.ZOO),
                      c("POM" = 1,"ALG" = param$f.POM),
                      c("DOM" = 1,"ALG" = param$f.DOM)) 

nu.gro.ZOO    <- calc.stoich.coef(alpha       = alpha,
                                  name        = "gro.ZOO",
                                  subst       = subst.gro.ZOO,
                                  subst.norm  = "ZOO",
                                  nu.norm     = 1,
                                  constraints = const.gro.ZOO)

# respiration of zooplankton:

subst.resp.ZOO <- c("NH4","HPO4","HCO3","O2","H","H2O","ZOO")

nu.resp.ZOO    <- calc.stoich.coef(alpha       = alpha,
                                   name        = "resp.ZOO",
                                   subst       = subst.resp.ZOO,
                                   subst.norm  = "ZOO",
                                   nu.norm     = -1)

# death of zooplankton:

subst.death.ZOO <- c("NH4","HPO4","HCO3","O2","H","H2O","ZOO","POM")

const.death.ZOO <- list(c("ZOO" = min(param$alpha.N.ZOO/param$alpha.N.POM,
                                      param$alpha.P.ZOO/param$alpha.P.POM,
                                      1),
                          "POM" = 1))
               # this constraint guarantees that no P and N is required for
               # death; if the P or N content of POM is larger than of ZOO
               # part of the body mass is mineralized
               # care must still be taken regarding potential oxygen consum
               
nu.death.ZOO    <- calc.stoich.coef(alpha       = alpha,
                                    name        = "death.ZOO",
                                    subst       = subst.death.ZOO,
                                    subst.norm  = "ZOO",
                                    nu.norm     = -1,
                                    constraints = const.death.ZOO)

# growth of phytoplankton with ammonium as nitrogen source:

subst.gro.ALG.NH4 <- c("NH4","HPO4","HCO3","O2","H","H2O","ALG")

nu.gro.ALG.NH4    <- calc.stoich.coef(alpha       = alpha,
                                      name        = "gro.ALG.NH4",
                                      subst       = subst.gro.ALG.NH4,
                                      subst.norm  = "ALG",
                                      nu.norm     = 1)

# growth of phytoplankton with nitrate as nitrogen source:

subst.gro.ALG.NO3 <- c("NO3","HPO4","HCO3","O2","H","H2O","ALG")

nu.gro.ALG.NO3    <- calc.stoich.coef(alpha       = alpha,
                                      name        = "gro.ALG.NO3",
                                      subst       = subst.gro.ALG.NO3,
                                      subst.norm  = "ALG",
                                      nu.norm     = 1)

# respiration of phytoplankton:

subst.resp.ALG <- c("NH4","HPO4","HCO3","O2","H","H2O","ALG")

nu.resp.ALG    <- calc.stoich.coef(alpha       = alpha,
                                   name        = "resp.ALG",
                                   subst       = subst.resp.ALG,
                                   subst.norm  = "ALG",
                                   nu.norm     = -1)

# death of phytoplankton:

subst.death.ALG <- c("NH4","HPO4","HCO3","O2","H","H2O","ALG","POM")

const.death.ALG <- list(c("ALG" = min(param$alpha.N.ALG/param$alpha.N.POM,
                                      param$alpha.P.ALG/param$alpha.P.POM,
                                      1),
                          "POM" = 1))
               # this constraint guarantees that no P and N is required for
               # death; if the P or N content of POM is larger than of ALG
               # part of the body mass is mineralized
               # care must still be taken regarding potential oxygen consumption

nu.death.ALG    <- calc.stoich.coef(alpha       = alpha,
                                    name        = "death.ALG",
                                    subst       = subst.death.ALG,
                                    subst.norm  = "ALG",
                                    nu.norm     = -1,
                                    constraints = const.death.ALG)
                                 
# hydrolysis:

subst.hydrolysis <- c("NH4","HPO4","HCO3","O2","H","H2O","POM","DOM")

const.hydrol     <- list(c("POM" = min(param$alpha.N.POM/param$alpha.N.DOM,
                                       param$alpha.P.POM/param$alpha.P.DOM,
                                       1),
                           "DOM" = 1))
               # this constraint guarantees that no P and N is required for
               # hydrolysis; if the P or N content of DOM is larger than of POM
               # part of the organic particle is mineralized

nu.hydrol <- calc.stoich.coef(alpha       = alpha,
                              name        = "hydrolysis",
                              subst       = subst.hydrolysis,
                              subst.norm  = "POM",
                              nu.norm     = -1,
                              constraints = const.hydrol)
                                 
# oxic mineralization of dissolved organic matter:

subst.miner.ox <- c("NH4","HPO4","HCO3","O2","H","H2O","DOM")

nu.miner.ox    <- calc.stoich.coef(alpha       = alpha,
                                   name        = "miner.ox",
                                   subst       = subst.miner.ox,
                                   subst.norm  = "DOM",
                                   nu.norm     = -1)

# nitrification:

subst.nitri <- c("NH4","NO3","O2","H","H2O")

nu.nitri    <- calc.stoich.coef(alpha       = alpha,
                                name        = "nitri",
                                subst       = subst.nitri,
                                subst.norm  = "NH4",
                                nu.norm     = -1)


# compile and print stoichiometric matrix:
# ----------------------------------------

nu <- rbind(nu.gro.ALG.NH4,
            nu.gro.ALG.NO3,
            nu.resp.ALG,
            nu.death.ALG,
            nu.gro.ZOO,
            nu.resp.ZOO,
            nu.death.ZOO,
            nu.hydrol,
            nu.miner.ox,
            nu.nitri)

print(nu,digits=2)


# write composition matrix and stoichiometric matrix to files:
# ============================================================

write.table(alpha,file="stoichcalc_11_exercise2_alpha.dat",sep="\t",col.names=NA)
write.table(nu,file="stoichcalc_11_exercise2_nu.dat",sep="\t",col.names=NA)

