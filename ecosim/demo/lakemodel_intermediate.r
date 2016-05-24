# ==============================================================================
# Lake Model of Intermediate Complexity
# ==============================================================================


# Load package:
# =============

library(ecosim)


# Definition of parameters:
# =========================

param <- list(alpha.O.ALG     = 0.50,       # gO/gALG
            alpha.H.ALG     = 0.07,         # gH/gALG
            alpha.N.ALG     = 0.06,         # gN/gALG
            alpha.P.ALG     = 0.005,        # gP/gALG
            alpha.O.ZOO     = 0.50,         # gO/gZOO
            alpha.H.ZOO     = 0.07,         # gH/gZOO
            alpha.N.ZOO     = 0.06,         # gN/gZOO
            alpha.P.ZOO     = 0.01,         # gP/gZOO
            alpha.O.POM     = 0.39,         # gO/gPOM
            alpha.H.POM     = 0.07,         # gH/gPOM
            alpha.N.POM     = 0.06,         # gN/gPOM
            alpha.P.POM     = 0.007,        # gP/gPOM
            Y.ZOO           = 0.2,          # gZOO/gALG
            f.e             = 0.2,          # gPOM/gALG
            k.gro.ALG       = 0.9,          # 1/d 
            k.gro.ZOO       = 0.6,          # m3/gDM/d
            k.resp.ALG      = 0.10,         # 1/d
            k.resp.ZOO      = 0.10,         # 1/d
            k.death.ALG     = 0.10,         # 1/d
            k.death.ZOO     = 0.08,         # 1/d
            k.miner.POM     = 0.02,         # 1/d
            k.miner.POM.sed = 5.0,          # gDM/m2/d
            K.POM.miner.sed = 10,           # gDM/m2
            K.HPO4          = 0.002,        # gP/m3
            K.O2.miner      = 0.5,          # gO/m3
            A               = 8.5e+006,     # m2
            h.epi           = 4,            # m
            h.hypo          = 8,            # m
            Q.in            = 4,            # m3/s
            C.ALG.ini       = 0.05,         # gDM/m3
            C.ZOO.ini       = 0.1,          # gDM/m3
            C.POM.ini       = 0.0,          # gDM/m3
            C.HPO4.ini      = 0.02,         # gP/m3
            C.O2.ini        = 10,           # gO/m3
            C.HPO4.in       = 0.04,         # gP/m3
            C.O2.in         = 10,           # gO/m3
            D.POM.ini       = 20,           # gDM/m2
            beta.ALG        = 0.046,        # 1/degC
            beta.ZOO        = 0.08,         # 1/degC
            beta.BAC        = 0.046,        # 1/degC
            T0              = 20,           # degC
            K.I             = 30,           # W/m2
            lambda.1        = 0.05,         # 1/m
            lambda.2        = 0.03,         # m2/gDM  
            v.ex.O2         = 1,            # m/d
            v.sed.POM       = 1,            # m/d 
            Kz.summer       = 0.02,         # m2/d
            Kz.winter       = 20,           # m2/d
            h.meta          = 5,            # m
            t.max           = 230,          # d
            I0.min          = 25,           # W/m2
            I0.max          = 225,          # W/m2
            T.min           = 5,            # degC
            T.max           = 25,           # degC
            p               = 101325)       # Pa

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
                   C.HPO4   = HPO4,
                   C.HCO3   = HCO3,
                   C.O2     = O2,
                   C.H      = H,
                   C.H2O    = H2O,
                   C.ALG    = ALG,
                   C.ZOO    = ZOO,
                   C.POM    = POM,
                   D.POM    = POM)

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
                                    "C.H","C.H2O","C.ALG","C.POM"),
                    subst.norm  = "C.ALG",
                    nu.norm     = -1,
                    constraints = list(c("C.ALG" = param$Y.ALG.death,
                                         "C.POM" = 1)))

# Growth of zooplankton:
# ----------------------

nu.gro.ZOO <- 
   calc.stoich.coef(alpha       = alpha,
                    name        = "gro.ZOO",
                    subst       = c("C.NH4","C.HPO4","C.HCO3","C.O2",
                                    "C.H","C.H2O","C.ALG","C.ZOO","C.POM"),
                    subst.norm  = "C.ZOO",
                    nu.norm     = 1,
                    constraints = list(c("C.ZOO" = 1,
                                         "C.ALG" = param$Y.ZOO),
                                       c("C.POM" = 1,
                                         "C.ALG" = param$f.e)))

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
                                    "C.H","C.H2O","C.ZOO","C.POM"),
                    subst.norm  = "C.ZOO",
                    nu.norm     = -1,
                    constraints = list(c("C.ZOO" = param$Y.ZOO.death,
                                         "C.POM" = 1)))

# Mineralization of suspended organic particles:
# ----------------------------------------------

nu.miner.POM <- 
   calc.stoich.coef(alpha       = alpha,
                    name        = "miner.POM",
                    subst       = c("C.NH4","C.HPO4","C.HCO3","C.O2",
                                    "C.H","C.H2O","C.POM"),
                    subst.norm  = "C.POM",
                    nu.norm     = -1)

# Mineralization of sedimented organic particles:
# -----------------------------------------------

nu.miner.POM.sed <- 
   calc.stoich.coef(alpha       = alpha,
                    name        = "miner.POM.sed",
                    subst       = c("C.NH4","C.HPO4","C.HCO3","C.O2",
                                    "C.H","C.H2O","D.POM"),
                    subst.norm  = "D.POM",
                    nu.norm     = -1)

# Sedimentation organic particles:
# --------------------------------

nu.sed.POM <- 
   calc.stoich.coef(alpha       = alpha,
                    name        = "sed.POM",
                    subst       = c("C.POM","D.POM"),
                    subst.norm  = "C.POM",
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
            nu.miner.POM,
            nu.miner.POM.sed,
            nu.sed.POM)
print(nu)


# Definition of transformation processes:
# =======================================

gro.ALG <- 
   new(Class  = "process",
       name   = "gro.ALG",
       rate   = expression(k.gro.ALG
                           *exp(beta.ALG*(T-T0))
                           *(C.HPO4/(K.HPO4+C.HPO4))
                           *log((K.I+I0)/
                                (K.I+I0*exp(-(lambda.1+lambda.2*C.ALG)*h.epi)))
                            /((lambda.1+lambda.2*C.ALG)*h.epi)
                           *C.ALG),
       stoich = as.list(nu["gro.ALG.NO3",]))

# Respiration of algae:
# ---------------------

resp.ALG <- 
   new(Class  = "process",
       name   = "resp.ALG",
       rate   = expression(k.resp.ALG*exp(beta.ALG*(T-T0))*C.ALG),
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
                           *C.ALG
                           *C.ZOO),
       stoich = as.list(nu["gro.ZOO",]))

# Respiration of zooplankton:
# ---------------------------

resp.ZOO <- 
   new(Class  = "process",
       name   = "resp.ZOO",
       rate   = expression(k.resp.ZOO*exp(beta.ZOO*(T-T0))*C.ZOO),
       stoich = as.list(nu["resp.ZOO",]))

# Death of zooplankton:
# ---------------------

death.ZOO <- 
   new(Class = "process",
       name   = "death.ZOO",
       rate   = expression(k.death.ZOO*C.ZOO),
       stoich = as.list(nu["death.ZOO",]))

# Mineralization of suspended organic particles:
# ----------------------------------------------

miner.POM <- 
   new(Class  = "process",
       name   = "miner.POM",
       rate   = expression(k.miner.POM
                           *exp(beta.BAC*(T-T0))
                           *C.O2/(K.O2.miner+C.O2)
                           *C.POM),
       stoich = as.list(nu["miner.POM",]))

# Mineralization of sedimented organic particles:
# -----------------------------------------------

miner.POM.sed <- 
   new(Class  = "process",
       name   = "miner.POM.sed",
       rate   = expression(k.miner.POM.sed
                           *exp(beta.BAC*(T-T0))
                           *C.O2/(K.O2.miner+C.O2)
                           *D.POM/(K.POM.miner.sed+D.POM)),
       stoich = as.list(nu["miner.POM.sed",]),
       pervol = F)

# Sedimentation:
# --------------

sed.POM <- 
   new(Class  = "process",
       name   = "sed.POM",
       rate   = expression(v.sed.POM/h.hypo*C.POM),
       stoich = as.list(nu["sed.POM",]))


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
Kz[i]	   <-eval(cond.gen$Kz, 	     envir=c(param, t=t[i]))
I0[i]	   <-eval(cond.epi$I0, 	     envir=c(param, t=t[i]))
T[i] 	   <-eval(cond.epi$T,  	     envir=c(param, t=t[i]))
C.O2.sat[i]<-eval(cond.epi$C.O2.sat, envir=c(param, T=T[i]))
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
                               C.O2   = expression(C.O2.ini),    # gN/m3
                               C.ALG  = expression(C.ALG.ini),   # gDM/m3
                               C.ZOO  = expression(C.ZOO.ini),   # gDM/m3
                               C.POM  = expression(C.POM.ini)),  # gDM/m3
       input            = list(C.O2   = expression(v.ex.O2*A
                                                   *(C.O2.sat-C.O2))), # gas ex.
       inflow           = expression(Q.in*86400),                # m3/d
       inflow.conc      = list(C.HPO4 = expression(C.HPO4.in),
                               C.O2   = expression(C.O2.in),
                               C.ALG  = 0,
                               C.ZOO  = 0,
                               C.POM  = 0),
       outflow          = expression(Q.in*86400),
       cond             = cond.epi,
       processes        = list(gro.ALG,resp.ALG,death.ALG,
                               gro.ZOO,resp.ZOO,death.ZOO,
                               miner.POM))

# Hypolimnion:
# ------------

hypolimnion <- 
   new(Class            = "reactor",
       name             = "Hypo",
       volume.ini       = expression(A*h.hypo),
       area             = expression(A),
       conc.pervol.ini  = list(C.HPO4 = expression(C.HPO4.ini),  # gP/m3
                               C.O2   = expression(C.O2.ini),    # gN/m3
                               C.ALG  = expression(C.ALG.ini),   # gDM/m3
                               C.ZOO  = expression(C.ZOO.ini),   # gDM/m3
                               C.POM  = expression(C.POM.ini)),  # gDM/m3
       conc.perarea.ini = list(D.POM  = expression(D.POM.ini)),  # gDM/m2
       cond             = cond.hypo,
       processes        = list(gro.ALG,resp.ALG,death.ALG,
                               gro.ZOO,resp.ZOO,death.ZOO,
                               miner.POM,miner.POM.sed,sed.POM))


# Definition of links:
# ====================

# Exchange between epilimnion and hypolimnion:
# --------------------------------------------

metalimnion <-
   new(Class      = "link",
       name       = "Metalimnion",
       from       = "Epi",
       to         = "Hypo",
       qadv.spec  = list(C.POM = expression(v.sed.POM*A)),
       qdiff.gen  = expression(A/h.meta*Kz))


# Definition of system:
# =====================

# Lake system:
# ------------

system <- new(Class    = "system",
              name     = "Lake",
              reactors = list(epilimnion,hypolimnion),
              links    = list(metalimnion),
              cond     = cond.gen,
              param    = param,
              t.out    = seq(0,365,by=1))    #number of days  365, 730


# Perform simulation:
# ===================

res <- calcres(system)


# Plot results:
# =============
                 
#plot.res(res)              # plot to screen

#plotres(res,colnames=list(c("C.HPO4.Epi","C.HPO4.Hypo"),                          
#                          c("C.O2.Epi","C.O2.Hypo"),                              
#                          c("C.ALG.Epi","C.ALG.Hypo"),                            
#                          c("C.ZOO.Epi","C.ZOO.Hypo"),                            
#                          c("C.POM.Epi","C.POM.Hypo"),                            
#                            "D.POM.Hypo"))                                          

plotres(res      = res,
        colnames = list(c("C.HPO4.Epi","C.HPO4.Hypo"),                          
                        c("C.O2.Epi","C.O2.Hypo"),                              
                        c("C.ALG.Epi","C.ALG.Hypo"),                            
                        c("C.ZOO.Epi","C.ZOO.Hypo"),                            
                        c("C.POM.Epi","C.POM.Hypo"),                            
                          "D.POM.Hypo"),
        file     = "lakemodel_intermediate.pdf",
        width    = 8,
        height   = 6)
                                                                                   
                           