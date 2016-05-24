# ==============================================================================
# STOICHCALC - R-Routines for Solving Stoichiometric Equations
# ==============================================================================
#
# Example corresponding to section 5.1 of the paper cited below       22.02.2010
# =============================================================
#
# Literature: Peter Reichert and Nele Schuwirth
#             A generic framework for deriving process stoichiometry
#             in environmental models
#             submitted to Environmental Modelling and Software
#
# ==============================================================================


# load library functions:
# =======================

library("stoichcalc")


# Define composition:
# ===================

# define composition parameters:
# ------------------------------

param <- list(a.N.ALG = 0.06, # gN/gALG
              a.P.ALG = 0.01) # gP/gALG


# define composition vectors of substances:
# -----------------------------------------

NH4   <- c(N = 1)             # gN/gNH4-N
HPO4  <- c(P = 1)             # gP/gHPO4-P
ALG   <- c(N = param$a.N.ALG, # gN/gALG
           P = param$a.P.ALG) # gP/gALG


# define list of composition vectors:
# -----------------------------------

subst.comp <- list(NH4  = NH4,
                   HPO4 = HPO4,
                   ALG  = ALG)


# compile and print composition matrix:
# -------------------------------------

alpha <- calc.comp.matrix(subst.comp)

print(alpha,digits=2)


# Derivation of Process Stoichiometry:
# ====================================

# calculate stoichiometries of growth and respiration processes:
# --------------------------------------------------------------

# growth of algae:

subst.gro.ALG.NH4 <- c("NH4","HPO4","ALG")

basis <- calc.stoich.basis(alpha = alpha,
                           subst = subst.gro.ALG.NH4)

nu.gro.ALG.NH4 <- calc.stoich.coef(alpha      = alpha,
                                   name       = "gro.ALG.NH4",
                                   subst      = subst.gro.ALG.NH4,
                                   subst.norm = "ALG",
                                   nu.norm    = 1)

# respiration of algae:

subst.resp.ALG <- c("NH4","HPO4","ALG")

basis <- calc.stoich.basis(alpha = alpha,
                           subst = subst.resp.ALG)

nu.resp.ALG <- calc.stoich.coef(alpha      = alpha,
                                name       = "resp.ALG",
                                subst      = subst.resp.ALG,
                                subst.norm = "ALG",
                                nu.norm    = -1)


# compile and print stoichiometric matrix of all processes:
# ---------------------------------------------------------

nu <- rbind(nu.gro.ALG.NH4,
            nu.resp.ALG)

print(nu,digits=2)


