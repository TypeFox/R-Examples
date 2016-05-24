# ==============================================================================
# Complex Lake Model
# ==============================================================================


# Load package:
# =============

library(ecosim)


# Definition of parameters:
# =========================

param <- list(alpha.O.ALG        = 0.50,         # gO/gALG
            alpha.H.ALG          = 0.07,         # gH/gALG
            alpha.N.ALG          = 0.06,         # gN/gALG
            alpha.P.ALG          = 0.005,        # gP/gALG
            alpha.O.ZOO          = 0.50,         # gO/gZOO
            alpha.H.ZOO          = 0.07,         # gH/gZOO
            alpha.N.ZOO          = 0.06,         # gN/gZOO
            alpha.P.ZOO          = 0.01,         # gP/gZOO
            alpha.O.POM          = 0.39,         # gO/gPOM
            alpha.H.POM          = 0.07,         # gH/gPOM
            alpha.N.POM          = 0.06,         # gN/gPOM
            alpha.P.POM          = 0.007,        # gP/gPOM
            Y.ZOO                = 0.2,          # gZOO/gALG
            f.e                  = 0.2,          # gPOM/gALG
            f.I                  = 0.2,          # gPOMI/(gPOMD+gPOMI)
            k.gro.ALG            = 0.9,          # 1/d
            k.gro.ZOO            = 0.6,          # m3/gDM/d
            k.resp.ALG           = 0.10,         # 1/d
            k.resp.ZOO           = 0.10,         # 1/d
            k.death.ALG          = 0.10,         # 1/d
            k.death.ZOO          = 0.08,         # 1/d
            k.nitri              = 0.1,          # gN/m3/d
            k.miner.ox.POM       = 0.02,         # 1/d
            k.miner.ox.POM.sed   = 5.0,          # gDM/m2/d
            k.miner.anox.POM.sed = 5.0,          # gDM/m2/d
            K.POM.miner.sed      = 10,           # gDM/m2
            K.HPO4               = 0.002,        # gP/m3
            K.N                  = 0.04,         # gN/m3
            p.NH4                = 5,            # - 
            K.O2.ZOO             = 0.2,          # gO/m3
            K.O2.resp            = 0.5,          # gO/m3
            K.O2.nitri           = 0.4,          # gO/m3
            K.O2.miner           = 0.5,          # gO/m3
            K.NO3.miner          = 0.1,          # gN/m3
            K.NH4.nitri          = 0.5,          # gN/m3
            A                    = 8.5e+006,     # m2
            h.epi                = 4,            # m
            h.hypo               = 8,            # m
            Q.in                 = 4,            # m3/s
            C.ALG.ini            = 0.05,         # gDM/m3
            C.ZOO.ini            = 0.1,          # gDM/m3
            C.POMD.ini           = 0.0,          # gDM/m3
            C.POMI.ini           = 0.0,          # gDM/m3
            C.HPO4.ini           = 0.05,         # gP/m3
            C.NH4.ini            = 0.1,          # gN/m3
            C.NO3.ini            = 0.5,          # gN/m3
            C.O2.ini             = 10,           # gO/m3
            C.HPO4.in            = 0.04,         # gP/m3
            C.NO3.in             = 0.5,          # gN/m3
            C.O2.in              = 10,           # gO/m3
            D.POMD.ini           = 0,            # gDM/m2
            D.POMI.ini           = 0,            # gDM/m2
            beta.ALG             = 0.046,        # 1/degC
            beta.ZOO             = 0.08,         # 1/degC
            beta.BAC             = 0.046,        # 1/degC
            T0                   = 20,           # degC
            K.I                  = 30,           # W/m2
            lambda.1             = 0.05,         # 1/m
            lambda.2             = 0.03,         # m2/gDM  
            v.ex.O2              = 1,            # m/d
            v.sed.POM            = 1,            # m/d
            Kz.summer            = 0.02,         # m2/d
            Kz.winter            = 20,           # m2/d
            h.meta               = 5,            # m
            t.max                = 230,          # d
            I0.min               = 25,           # W/m2
            I0.max               = 225,          # W/m2
            T.min                = 5,            # degC
            T.max                = 25,           # degC
            p                    = 101325)       # Pa

# choose carbon fractions to guarantee that the fractions sum to unity:
param$alpha.C.ALG <- 1 - (param$alpha.O.ALG + param$alpha.H.ALG + 
                        param$alpha.N.ALG + param$alpha.P.ALG)
param$alpha.C.ZOO <- 1 - (param$alpha.O.ZOO + param$alpha.H.ZOO + 
                        param$alpha.N.ZOO + param$alpha.P.ZOO)
param$alpha.C.POM <- 1 - (param$alpha.O.POM + param$alpha.H.POM + 
                        param$alpha.N.POM + param$alpha.P.POM)

# choose yield of death to guarantee that no nutrients are required
# (oxygen content of POM was reduced to avoid need of oxygen):
param$Y.ALG.death <-  min(1,
                        param$alpha.N.ALG/param$alpha.N.POM,
                        param$alpha.P.ALG/param$alpha.P.POM,
                        param$alpha.C.ALG/param$alpha.C.POM)
param$Y.ZOO.death <-  min(1,
                        param$alpha.N.ZOO/param$alpha.N.POM,
                        param$alpha.P.ZOO/param$alpha.P.POM,
                        param$alpha.C.ZOO/param$alpha.C.POM)

# Construction of composition matrix:
# -----------------------------------

NH4    <- c(H      = 4*1/14,              # gH/gNH4-N
            N      = 1,                   # gN/gNH4-N
            charge = 1/14)                # chargeunits/gNH4-N
NO3    <- c(O      = 3*16/14,             # gO/gNO3-N
            N      = 1,                   # gN/gNO3-N
            charge = -1/14)               # chargeunits/gNO3-N
N2     <- c(N      = 1)                   # gN/gN2-N
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
H2O    <- c(O      = 1*16,                # gO/molH2O
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

subst.comp <- list(C.NH4    = NH4,
                   C.NO3    = NO3,
                   C.N2     = N2,
                   C.HPO4   = HPO4,
                   C.HCO3   = HCO3,
                   C.O2     = O2,
                   C.H      = H,
                   C.H2O    = H2O,
                   C.ALG    = ALG,
                   C.ZOO    = ZOO,
                   C.POMD   = POM,
                   D.POMD   = POM,
                   C.POMI   = POM,
                   D.POMI   = POM)

alpha <- calc.comp.matrix(subst.comp)
print(alpha)


# Derivation of Process Stoichiometry:
# ====================================

# Growth of algae on ammonium:
# ----------------------------

nu.gro.ALG.NH4 <- 
   calc.stoich.coef(alpha       = alpha,
                    name        = "gro.ALG.NH4",
                    subst       = c("C.NH4","C.HPO4","C.HCO3","C.O2",
                                    "C.H","C.H2O","C.ALG"),
                    subst.norm  = "C.ALG",
                    nu.norm     = 1)

# Growth of algae on nitrate:
# ---------------------------

nu.gro.ALG.NO3 <- 
   calc.stoich.coef(alpha       = alpha,
                    name        = "gro.ALG.NO3",
                    subst       = c("C.NO3","C.HPO4","C.HCO3","C.O2",
                                    "C.H","C.H2O","C.ALG"),
                    subst.norm  = "C.ALG",
                    nu.norm     = 1)

# Respiration of algae:
# ---------------------

nu.resp.ALG <- 
   calc.stoich.coef(alpha       = alpha,
                    name        = "resp.ALG",
                    subst       = c("C.NH4","C.HPO4","C.HCO3","C.O2",
                                    "C.H","C.H2O","C.ALG"),
                    subst.norm  = "C.ALG",
                    nu.norm     = -1)

# Death of algae:
# ---------------

nu.death.ALG <- 
   calc.stoich.coef(alpha       = alpha,
                    name        = "death.ALG",
                    subst       = c("C.NH4","C.HPO4","C.HCO3","C.O2",
                                    "C.H","C.H2O","C.ALG","C.POMD","C.POMI"),
                    subst.norm  = "C.ALG",
                    nu.norm     = -1,
                    constraints = list(c("C.ALG"  = param$Y.ALG.death,
                                         "C.POMD" = 1,
                                         "C.POMI" = 1),
                                       c("C.POMD" = -param$f.I,
                                         "C.POMI" = 1-param$f.I)))

# Growth of zooplankton:
# ----------------------

nu.gro.ZOO <- 
   calc.stoich.coef(alpha       = alpha,
                    name        = "gro.ZOO",
                    subst       = c("C.NH4","C.HPO4","C.HCO3","C.O2","C.H",
                                    "C.H2O","C.ALG","C.ZOO","C.POMD","C.POMI"),
                    subst.norm  = "C.ZOO",
                    nu.norm     = 1,
                    constraints = list(c("C.ZOO"  = 1,
                                         "C.ALG"  = param$Y.ZOO),
                                       c("C.POMD" = 1,
                                         "C.POMI" = 1,
                                         "C.ALG"  = param$f.e),
                                       c("C.POMD" = -param$f.I,
                                         "C.POMI" = 1-param$f.I)))

# Respiration of zooplankton:
# ---------------------------

nu.resp.ZOO <- 
   calc.stoich.coef(alpha       = alpha,
                    name        = "resp.ZOO",
                    subst       = c("C.NH4","C.HPO4","C.HCO3","C.O2",
                                    "C.H","C.H2O","C.ZOO"),
                    subst.norm  = "C.ZOO",
                    nu.norm     = -1)

# Death of zooplankton:
# ---------------------

nu.death.ZOO <- 
   calc.stoich.coef(alpha       = alpha,
                    name        = "death.ZOO",
                    subst       = c("C.NH4","C.HPO4","C.HCO3","C.O2",
                                    "C.H","C.H2O","C.ZOO","C.POMD","C.POMI"),
                    subst.norm  = "C.ZOO",
                    nu.norm     = -1,
                    constraints = list(c("C.ZOO"  = param$Y.ZOO.death,
                                         "C.POMD" = 1,
                                         "C.POMI" = 1),
                                       c("C.POMD" = -param$f.I,
                                         "C.POMI" = 1-param$f.I)))

# Nitrification:
# --------------

nu.nitri <- 
   calc.stoich.coef(alpha       = alpha,
                    name        = "nitri",
                    subst       = c("C.NH4","C.NO3","C.O2","C.H","C.H2O"),
                    subst.norm  = "C.NH4",
                    nu.norm     = -1)

# Oxic mineralization of suspended organic particles:
# ---------------------------------------------------

nu.miner.ox.POM <- 
   calc.stoich.coef(alpha       = alpha,
                    name        = "miner.ox.POM",
                    subst       = c("C.NH4","C.HPO4","C.HCO3","C.O2",
                                    "C.H","C.H2O","C.POMD"),
                    subst.norm  = "C.POMD",
                    nu.norm     = -1)

# Oxic mineralization of sedimented organic particles:
# ----------------------------------------------------

nu.miner.ox.POM.sed <- 
   calc.stoich.coef(alpha       = alpha,
                    name        = "miner.ox.POM.sed",
                    subst       = c("C.NH4","C.HPO4","C.HCO3","C.O2",
                                    "C.H","C.H2O","D.POMD"),
                    subst.norm  = "D.POMD",
                    nu.norm     = -1)

# Anoxic mineralization of sedimented organic particles:
# ------------------------------------------------------

nu.miner.anox.POM.sed <- 
   calc.stoich.coef(alpha       = alpha,
                    name        = "miner.anox.POM.sed",
                    subst       = c("C.NH4","C.NO3","C.HPO4","C.HCO3","C.N2",
                                    "C.H","C.H2O","D.POMD"),
                    subst.norm  = "D.POMD",
                    nu.norm     = -1,
                    constraints = list(c("C.NO3" = 1,
                                         "C.N2"  = 1)))

# Sedimentation of degradable organic particles:
# ----------------------------------------------

nu.sed.POMD <- 
   calc.stoich.coef(alpha       = alpha,
                    name        = "sed.POMD",
                    subst       = c("C.POMD","D.POMD"),
                    subst.norm  = "C.POMD",
                    nu.norm     = -1)

# Sedimentation of inert organic particles:
# -----------------------------------------

nu.sed.POMI <- 
   calc.stoich.coef(alpha       = alpha,
                    name        = "sed.POMI",
                    subst       = c("C.POMI","D.POMI"),
                    subst.norm  = "C.POMI",
                    nu.norm     = -1)

# Combine process stoichiometries to stoichiometric matrix:
# ---------------------------------------------------------

nu <- rbind(nu.gro.ALG.NH4,
            nu.gro.ALG.NO3,
            nu.resp.ALG,
            nu.death.ALG,
            nu.gro.ZOO,
            nu.resp.ZOO,
            nu.death.ZOO,
            nu.nitri,
            nu.miner.ox.POM,
            nu.miner.ox.POM.sed,
            nu.miner.anox.POM.sed,
            nu.sed.POMD,
            nu.sed.POMI)
print(nu)


# Definition of transformation processes:
# =======================================

# Growth of algae on ammonium:
# ----------------------------

gro.ALG.NH4 <- 
   new(Class  = "process",
       name   = "gro.ALG.NH4",
       rate   = expression(k.gro.ALG
                            *exp(beta.ALG*(T-T0))
                            *log((K.I+I0)/
                                (K.I+I0*exp(-(lambda.1+lambda.2*C.ALG)*h.epi)))
                                /((lambda.1+lambda.2*C.ALG)*h.epi)
                            *min(C.HPO4/(K.HPO4+C.HPO4),
                                (C.NH4+C.NO3)/(K.N+C.NH4+C.NO3))
                            *(p.NH4*C.NH4/(p.NH4*C.NH4+C.NO3))
                            *C.ALG),
       stoich = as.list(nu["gro.ALG.NH4",]))

# Growth of algae on nitrate:
# ---------------------------

gro.ALG.NO3 <- 
   new(Class  = "process",
       name   = "gro.ALG.NO3",
       rate   = expression(k.gro.ALG
                           *exp(beta.ALG*(T-T0))
                           *log((K.I+I0)/
                                (K.I+I0*exp(-(lambda.1+lambda.2*C.ALG)*h.epi)))
                                /((lambda.1+lambda.2*C.ALG)*h.epi)
                           *min(C.HPO4/(K.HPO4+C.HPO4),
                                (C.NH4+C.NO3)/(K.N+C.NH4+C.NO3))
                           *(C.NO3/(p.NH4*C.NH4+C.NO3))
                           *C.ALG),
       stoich = as.list(nu["gro.ALG.NO3",]))
       
# Respiration of algae:
# ---------------------

resp.ALG <- 
   new(Class  = "process",
       name   = "resp.ALG",
       rate   = expression(k.resp.ALG*exp(beta.ALG*(T-T0))*(C.O2/(K.O2.resp+C.O2))*C.ALG),
       stoich = as.list(nu["resp.ALG",]))

# Death of algae:
# ---------------

death.ALG <- 
   new(Class = "process",
       name   = "death.ALG",
       rate   = expression(k.death.ALG*C.ALG),
       stoich = as.list(nu["death.ALG",]))

# Growth of zooplankton:
# ----------------------

gro.ZOO <- 
   new(Class  = "process",
       name   = "gro.ZOO",
       rate   = expression(k.gro.ZOO*exp(beta.ZOO*(T-T0))
                           *(C.O2/(K.O2.ZOO+C.O2))
                           *C.ALG
                           *C.ZOO),
       stoich = as.list(nu["gro.ZOO",]))

# Respiration of zooplankton:
# ---------------------------

resp.ZOO <- 
   new(Class  = "process",
       name   = "resp.ZOO",
       rate   = expression(k.resp.ZOO*exp(beta.ZOO*(T-T0))*(C.O2/(K.O2.resp+C.O2))*C.ZOO),
       stoich = as.list(nu["resp.ZOO",]))

# Death of zooplankton:
# ---------------------

death.ZOO <- 
   new(Class = "process",
       name   = "death.ZOO",
       rate   = expression(k.death.ZOO*C.ZOO),
       stoich = as.list(nu["death.ZOO",]))

# Nitrification:
# --------------

nitri <- 
   new(Class = "process",
       name   = "nitri",
       rate   = expression(k.nitri
                           *exp(beta.BAC*(T-T0))
                           *min(C.NH4/(K.NH4.nitri+C.NH4),C.O2/(K.O2.nitri+C.O2))),
       stoich = as.list(nu["nitri",]))

# Oxic mineralization of suspended organic particles:
# ---------------------------------------------------

miner.ox.POM <- 
   new(Class  = "process",
       name   = "miner.ox.POM",
       rate   = expression(k.miner.ox.POM
                           *exp(beta.BAC*(T-T0))
                           *C.O2/(K.O2.miner+C.O2)
                           *C.POMD),
       stoich = as.list(nu["miner.ox.POM",]))

# Oxic mineralization of sedimented organic particles:
# ----------------------------------------------------

miner.ox.POM.sed <- 
   new(Class  = "process",
       name   = "miner.ox.POM.sed",
       rate   = expression(k.miner.ox.POM.sed
                           *exp(beta.BAC*(T-T0))
                           *C.O2/(K.O2.miner+C.O2)
                           *D.POMD/(K.POM.miner.sed+D.POMD)),
       stoich = as.list(nu["miner.ox.POM.sed",]),
       pervol = F)

# Anoxic mineralization of sedimented organic particles:
# ------------------------------------------------------

miner.anox.POM.sed <- 
   new(Class  = "process",
       name   = "miner.anox.POM.sed",
       rate   = expression(k.miner.anox.POM.sed
                           *exp(beta.BAC*(T-T0))
                           *C.NO3/(K.NO3.miner+C.NO3)
                           *(D.POMD/(K.POM.miner.sed+D.POMD))^2),
       stoich = as.list(nu["miner.anox.POM.sed",]),
       pervol = F)

# Sedimentation of degradable organic particles:
# ----------------------------------------------

sed.POMD <- 
   new(Class  = "process",
       name   = "sed.POMD",
       rate   = expression(v.sed.POM/h.hypo*C.POMD),
       stoich = as.list(nu["sed.POMD",]))

# Sedimentation of inert organic particles:
# -----------------------------------------

sed.POMI <- 
   new(Class  = "process",
       name   = "sed.POMI",
       rate   = expression(v.sed.POM/h.hypo*C.POMI),
       stoich = as.list(nu["sed.POMI",]))


# Definition of environmental conditions:
# =======================================

cond.epi  <- list(I0       = expression( 0.5*(I0.min+I0.max)
                                        +0.5*(I0.max-I0.min)
                                            *cos(2*pi/365.25*(t-t.max))),
                  T        = expression( 0.5*(T.min+T.max)
                                        +0.5*(T.max-T.min)
                                            *cos(2*pi/365.25*(t-t.max))),
                  C.O2.sat = expression(exp(7.7117-1.31403*log(T+45.93))
                                        *p/101325))
                                          
cond.hypo <- list(I0 = 0,
                  T  = 5)

cond.gen  <- list(Kz=expression( 0.5*(Kz.summer+Kz.winter)
                                -0.5*(Kz.winter-Kz.summer)
                                    *sign(cos(2*pi/365.25*(t-t.max))+0.4)))
                                    
# Plot the environmental conditions
t <-  1:365 
  
Kz <- numeric(0)
I0 <- numeric(0)
T <- numeric(0)  
C.O2.sat <- numeric(0)      
                       
for(i in 1:length(t))
{
   Kz[i]       <- eval(cond.gen$Kz, envir=c(param, t=t[i]))
   I0[i]       <- eval(cond.epi$I0, envir=c(param, t=t[i]))
   T[i]        <- eval(cond.epi$T, envir=c(param, t=t[i]))
   C.O2.sat[i] <- eval(cond.epi$C.O2.sat, envir=c(param, T=T[i]))
}

#Plots
par(mfrow=c(2,2),xaxs="i",yaxs="i",mar=c(4.5,4.5,2,1.5)+0.1)                                 
plot(t, Kz , ylim=c(0,25))
plot(t, I0, type="l")
plot(t, T, type="l")
plot(t, C.O2.sat, type="l")


# Definition of reactors:
# =======================

# Epilimnion:
# -----------

epilimnion <- 
   new(Class            = "reactor",
       name             = "Epi",
       volume.ini       = expression(A*h.epi),
       conc.pervol.ini  = list(C.HPO4 = expression(C.HPO4.ini),  # gP/m3
                               C.NH4  = expression(C.NH4.ini),   # gN/m3
                               C.NO3  = expression(C.NO3.ini),   # gN/m3
                               C.O2   = expression(C.O2.ini),    # gO/m3
                               C.ALG  = expression(C.ALG.ini),   # gDM/m3
                               C.ZOO  = expression(C.ZOO.ini),   # gDM/m3
                               C.POMD = expression(C.POMD.ini),  # gDM/m3
                               C.POMI = expression(C.POMI.ini)), # gDM/m3
       input            = list(C.O2   = expression(v.ex.O2*A
                                                   *(C.O2.sat-C.O2))), # gas ex.
       inflow           = expression(Q.in*86400),                # m3/d
       inflow.conc      = list(C.HPO4 = expression(C.HPO4.in),
                               C.NH4  = 0,
                               C.NO3  = expression(C.NO3.in),
                               C.O2   = expression(C.O2.in),
                               C.ALG  = 0,
                               C.ZOO  = 0,
                               C.POMD = 0,
                               C.POMI = 0),
       outflow          = expression(Q.in*86400),
       cond             = cond.epi,
       processes        = list(gro.ALG.NH4,gro.ALG.NO3,resp.ALG,death.ALG,
                               gro.ZOO,resp.ZOO,death.ZOO,
                               nitri,miner.ox.POM))

# Hypolimnion:
# ------------

hypolimnion <- 
   new(Class            = "reactor",
       name             = "Hypo",
       volume.ini       = expression(A*h.hypo),
       area             = expression(A),
       conc.pervol.ini  = list(C.HPO4 = expression(C.HPO4.ini),  # gP/m3
                               C.NH4  = expression(C.NH4.ini),   # gN/m3
                               C.NO3  = expression(C.NO3.ini),   # gN/m3
                               C.O2   = expression(C.O2.ini),    # gO/m3
                               C.ALG  = expression(C.ALG.ini),   # gDM/m3
                               C.ZOO  = expression(C.ZOO.ini),   # gDM/m3
                               C.POMD = expression(C.POMD.ini),  # gDM/m3
                               C.POMI = expression(C.POMI.ini)), # gDM/m3
       conc.perarea.ini = list(D.POMD = expression(D.POMD.ini),  # gDM/m2
                               D.POMI = expression(D.POMI.ini)), # gDM/m2
       cond             = cond.hypo,
       processes        = list(resp.ALG,death.ALG,
                               gro.ZOO,resp.ZOO,death.ZOO,
                               nitri,
                               miner.ox.POM,miner.ox.POM.sed,miner.anox.POM.sed,
                               sed.POMD,sed.POMI))


# Definition of links:
# ====================

# Exchange between epilimnion and hypolimnion:
# --------------------------------------------

metalimnion <-
   new(Class      = "link",
       name       = "Metalimnion",
       from       = "Epi",
       to         = "Hypo",
       qadv.spec  = list(C.POMD = expression(v.sed.POM*A),
                         C.POMI = expression(v.sed.POM*A)),
       qdiff.gen  = expression(A/h.meta*Kz))


#Definition of system:
#=====================

# Lake system:
# ------------

system <- new(Class    = "system",
              name     = "Lake",
              reactors = list(epilimnion,hypolimnion),
              links    = list(metalimnion),
              cond     = cond.gen,
              param    = param,
              t.out    = seq(0,730,by=1)) #number of days 365, 730, 1095 ,1460


# Perform simulation:
# ===================

res <- calcres(system)


# Plot results:
# =============
                 
#plotres(res)              # plot to screen

#plotres(res,colnames=
#             list(c("C.NH4.Epi","C.NH4.Hypo"),
#                  c("C.NO3.Epi","C.NO3.Hypo"),
#                  c("C.HPO4.Epi","C.HPO4.Hypo"),
#                  c("C.O2.Epi","C.O2.Hypo"),
#                  c("C.ALG.Epi","C.ALG.Hypo"),
#                  c("C.ZOO.Epi","C.ZOO.Hypo"),
#                  c("C.POMD.Epi","C.POMD.Hypo","C.POMI.Epi","C.POMI.Hypo"),
#                  "D.POMD.Hypo",
#                  "D.POMI.Hypo"))
                  
plotres(res,colnames=
             list(c("C.NH4.Epi","C.NH4.Hypo"),
                  c("C.NO3.Epi","C.NO3.Hypo"),
                  c("C.HPO4.Epi","C.HPO4.Hypo"),
                  c("C.O2.Epi","C.O2.Hypo"),
                  c("C.ALG.Epi","C.ALG.Hypo"),
                  c("C.ZOO.Epi","C.ZOO.Hypo"),
                  c("C.POMD.Epi","C.POMD.Hypo","C.POMI.Epi","C.POMI.Hypo"),
                  "D.POMD.Hypo",
                  "D.POMI.Hypo"),
         title    = "Complex Lake Model",
         file     = "lakemodel_complex.pdf",
         width    = 8,
         height   = 6)


