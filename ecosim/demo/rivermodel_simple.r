# ==============================================================================
# Simple River Model
# ==============================================================================


# Load package:
# =============

library(ecosim)


# Definition of parameters:
# =========================

param <- list(

#stoichiometric parameters
            alpha.O.ALG     = 0.50,         # gO/gALG
            alpha.H.ALG     = 0.07,         # gH/gALG
            alpha.N.ALG     = 0.06,         # gN/gALG
            alpha.P.ALG     = 0.01,         # gP/gALG
            alpha.O.POM     = 0.39,         # gO/gPOM
            alpha.H.POM     = 0.07,         # gH/gPOM
            alpha.N.POM     = 0.04,         # gN/gPOM
            alpha.P.POM     = 0.01,         # gP/gPOM

#kinetic parameters           
            k.gro.ALG       = 1.5,          # 1/d
            k.resp.ALG      = 0.10,         # 1/d
            k.nitri1        = 10.0,         # gNH4-N/m3/d
            k.nitri2        = 10.0,         # gNO2-N/m3/d
            k.miner.POM     = 0.5,          # 1/d
            K.HPO4.ALG      = 0.002,        # gP/m3
            K.N.ALG         = 0.04,         # gN/m3
            p.NH4           = 5,            # -
            K.NH4.nitri     = 0.5,          # gN/m3
            K.NO2.nitri     = 0.5,          # gN/m3
            K.O2.nitri      = 0.5,          # gO/m3
            K.O2.resp       = 0.5,          # gO/m3
            K.O2.miner      = 0.5,          # gO/m3
            beta.ALG        = 0.046,        # 1/degC
            beta.BAC        = 0.046,        # 1/degC
            beta.nitri1     = 0.098,        # 1/degC
            beta.nitri2     = 0.069,        # 1/degC
            K.I             = 200,          # W/m2
            lambda          = 0.1,          # 1/m
            K2.O2           = 10,           # 1/d

#River geometry        
            L               = 2000,         # m
            w               = 20,           # m
            h               = 0.5,          # m    
                   
#Input and initial conditions
            Q.in            = 4,            # m3/s
            D.ALG           = 50,           # gDM/m2
            D.POM           = 50,           # gDM/m2
            C.HPO4.ini      = 0.4,          # gP/m3
            C.NH4.ini       = 0.4,          # gN/m3
            C.NO3.ini       = 4,            # gN/m3
            C.NO2.ini       = 0,            # gN/m3
            C.O2.ini        = 10,           # gO/m3
            C.HPO4.in       = 0.4,          # gP/m3
            C.NH4.in        = 0.4,          # gN/m3
            C.NO3.in        = 4.0,          # gN/m3
            C.O2.in         = 10,           # gO/m3
            
#environmental conditions
            T0              = 20,           # degC
            T.min           = 15,           # degC
            T.max           = 25,           # degC
	      I0.max          = 600,          # W/m2
	      t.max.I         = 0.5,          # d
            t.max.T         = 0.6,          # d
	      p               = 101325)       # Pa

# choose carbon fractions to guarantee that the fractions sum to unity:
param$alpha.C.ALG <- 1 - (param$alpha.O.ALG + param$alpha.H.ALG + 
                          param$alpha.N.ALG + param$alpha.P.ALG)
param$alpha.C.POM <- 1 - (param$alpha.O.POM + param$alpha.H.POM + 
                          param$alpha.N.POM + param$alpha.P.POM)


# Construction of composition matrix:
# -----------------------------------

NH4    <- c(H      = 4*1/14,              # gH/gNH4-N      
            N      = 1,                   # gN/gNH4-N
            charge = 1/14)                # chargeunits/gNH4-N
NO2    <- c(O      = 2*16/14,             # gO/gNO2-N
            N      = 1,                   # gN/gNO2-N
            charge = -1/14)               # chargeunits/gNO2-N
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
POM    <- c(C      = param$alpha.C.POM,     # gC/gPOM
            O      = param$alpha.O.POM,     # gO/gPOM
            H      = param$alpha.H.POM,     # gH/gPOM
            N      = param$alpha.N.POM,     # gN/gPOM
            P      = param$alpha.P.POM)     # gP/gPOM

subst.comp <- list(C.NH4    = NH4,
                   C.NO2    = NO2,
                   C.NO3    = NO3,
                   C.HPO4   = HPO4,
                   C.HCO3   = HCO3,
                   C.O2     = O2,
                   C.H      = H,
                   C.H2O    = H2O,
                   D.ALG    = ALG,
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
                                    "C.H","C.H2O","D.ALG"),
                    subst.norm  = "D.ALG",
                    nu.norm     = 1)

# Growth of algae on nitrate:
# ---------------------------

nu.gro.ALG.NO3 <- 
   calc.stoich.coef(alpha       = alpha,
                    name        = "gro.ALG.NO3",
                    subst       = c("C.NO3","C.HPO4","C.HCO3","C.O2",
                                    "C.H","C.H2O","D.ALG"),
                    subst.norm  = "D.ALG",
                    nu.norm     = 1)

# Respiration of algae:
# ---------------------

nu.resp.ALG <- 
   calc.stoich.coef(alpha       = alpha,
                    name        = "resp.ALG",
                    subst       = c("C.NH4","C.HPO4","C.HCO3","C.O2",
                                    "C.H","C.H2O","D.ALG"),
                    subst.norm  = "D.ALG",
                    nu.norm     = -1)

# First step of nitrification:
# ----------------------------

nu.nitri1 <- 
   calc.stoich.coef(alpha       = alpha,
                    name        = "nitri1",
                    subst       = c("C.NH4","C.NO2","C.O2",
                                    "C.H","C.H2O"),
                    subst.norm  = "C.NH4",
                    nu.norm     = -1)

# Second step of nitrification:
# -----------------------------

nu.nitri2 <- 
   calc.stoich.coef(alpha       = alpha,
                    name        = "nitri2",
                    subst       = c("C.NO2","C.NO3","C.O2",
                                    "C.H","C.H2O"),
                    subst.norm  = "C.NO2",
                    nu.norm     = -1)

# Mineralization of organic particles:
# ------------------------------------

nu.miner.POM <- 
   calc.stoich.coef(alpha       = alpha,
                    name        = "miner.POM",
                    subst       = c("C.NH4","C.HPO4","C.HCO3","C.O2",
                                    "C.H","C.H2O","D.POM"),
                    subst.norm  = "D.POM",
                    nu.norm     = -1)

# Combine process stoichiometries to stoichiometric matrix:
# ---------------------------------------------------------

nu <- rbind(nu.gro.ALG.NH4,
            nu.gro.ALG.NO3,
            nu.resp.ALG,
            nu.nitri1,
            nu.nitri2,
            nu.miner.POM)
print(nu)

#write.table(nu,file="river1_nu.dat",sep="\t",col.names=NA)


# Definition of transformation processes:
# =======================================

# Growth of algae on ammonium:
# ----------------------------

gro.ALG.NH4 <- 
   new(Class  = "process",
       name   = "gro.ALG.NH4",
       rate   = expression(k.gro.ALG
                           *exp(beta.ALG*(T-T0))
                           *I0*exp(-lambda*h)/(K.I+I0*exp(-lambda*h))
                           *min(C.HPO4/(K.HPO4.ALG+C.HPO4),
                                (C.NH4+C.NO3)/(K.N.ALG+C.NH4+C.NO3))
                           *(p.NH4*C.NH4/(p.NH4*C.NH4+C.NO3))
                           *D.ALG),
       stoich = as.list(nu["gro.ALG.NH4",]),
       pervol = F)

# Growth of algae on nitrate:
# ---------------------------

gro.ALG.NO3 <- 
   new(Class  = "process",
       name   = "gro.ALG.NO3",
       rate   = expression(k.gro.ALG
                           *exp(beta.ALG*(T-T0))
                           *I0*exp(-lambda*h)/(K.I+I0*exp(-lambda*h))
                           *min(C.HPO4/(K.HPO4.ALG+C.HPO4),
                                (C.NH4+C.NO3)/(K.N.ALG+C.NH4+C.NO3))
                           *(C.NO3/(p.NH4*C.NH4+C.NO3))
                           *D.ALG),
       stoich = as.list(nu["gro.ALG.NO3",]),
       pervol = F)

# Respiration of algae:
# ---------------------

resp.ALG <- 
   new(Class  = "process",
       name   = "resp.ALG",
       rate   = expression(k.resp.ALG*exp(beta.ALG*(T-T0))
                           *C.O2/(K.O2.resp+C.O2)*D.ALG),
       stoich = as.list(nu["resp.ALG",]),
       pervol = F)

# First step of nitrification:
# ----------------------------

nitri1 <- 
   new(Class  = "process",
       name   = "nitri1",
       rate   = expression(k.nitri1
                           *exp(beta.nitri1*(T-T0))
                           *min(C.NH4/(K.NH4.nitri+C.NH4),
                                C.O2/(K.O2.nitri+C.O2))),
       stoich = as.list(nu["nitri1",]))

# Second step of nitrification:
# -----------------------------

nitri2 <- 
   new(Class  = "process",
       name   = "nitri2",
       rate   = expression(k.nitri2
                           *exp(beta.nitri2*(T-T0))
                           *min(C.NO2/(K.NO2.nitri+C.NO2),
                                C.O2/(K.O2.nitri+C.O2))),
       stoich = as.list(nu["nitri2",]))

# Mineralization of organic particles:
# ------------------------------------

miner.POM <- 
   new(Class  = "process",
       name   = "miner.POM",
       rate   = expression(k.miner.POM
                           *exp(beta.BAC*(T-T0))
                           *C.O2/(K.O2.miner+C.O2)
                           *D.POM),
       stoich = as.list(nu["miner.POM",]),
       pervol = F)


# Definition of environmental conditions:
# =======================================

cond <- list(I0       = expression(I0.max*0.5*(sign(cos(2*pi/1.0*(t-t.max.I)))+1)
                                   *cos(2*pi/1.0*(t-t.max.I))^2),
             T        = expression( 0.5*(T.min+T.max)
                                   +0.5*(T.max-T.min)
                                       *cos(2*pi/1.0*(t-t.max.T))),
             C.O2.sat = expression(exp(7.7117-1.31403*log(T+45.93))
                                   *p/101325))

# Plot the environmental conditions
t <-  seq(0,3,by=0.01) 

I0 <- numeric(0)
T <- numeric(0)  
C.O2.sat <- numeric(0) 

for(i in 1:length(t))
{
   I0[i]       <- eval(cond$I0, envir=c(param, t=t[i]))
   T[i]        <- eval(cond$T, envir=c(param, t=t[i]))
   C.O2.sat[i] <- eval(cond$C.O2.sat, envir=c(param, T=T[i]))
}

#Plots
par(mfrow=c(2,2),xaxs="i",yaxs="i",mar=c(4.5,4.5,2,1.5)+0.1)                                 
plot(t, I0, type="l")
plot(t, T, type="l")
plot(t, C.O2.sat, type="l")


# Definition of reactor:
# ======================

# River Sections:
# ---------------

reach1 <- 
   new(Class            = "reactor",
       name             = "R1",
       volume.ini       = expression(L*w*h),
       area             = expression(L*w),
       conc.pervol.ini  = list(C.HPO4 = expression(C.HPO4.ini),  # gP/m3
                               C.NH4  = expression(C.NH4.ini),   # gN/m3
                               C.NO2  = expression(C.NO2.ini),   # gN/m3
                               C.NO3  = expression(C.NO3.ini),   # gN/m3
                               C.O2   = expression(C.O2.ini)),   # gN/m3
       input            = list(C.O2   = expression(K2.O2*L*w*h
                                                   *(C.O2.sat-C.O2))), # gas exchange
       inflow           = expression(Q.in*86400),                # m3/d
       inflow.conc      = list(C.HPO4 = expression(C.HPO4.in),
                               C.NH4  = expression(C.NH4.in),
                               C.NO3  = expression(C.NO3.in),
                               C.O2   = expression(C.O2.sat)),
       processes        = list(gro.ALG.NH4,gro.ALG.NO3,resp.ALG, #     
                               nitri1,nitri2,miner.POM))

reach2 <- 
   new(Class            = "reactor",
       name             = "R2",
       volume.ini       = expression(L*w*h),
       area             = expression(L*w),
       conc.pervol.ini  = list(C.HPO4 = expression(C.HPO4.ini),  # gP/m3
                               C.NH4  = expression(C.NH4.ini),   # gN/m3
                               C.NO2  = expression(C.NO2.ini),   # gN/m3
                               C.NO3  = expression(C.NO3.ini),   # gN/m3
                               C.O2   = expression(C.O2.ini)),   # gN/m3
       input            = list(C.O2   = expression(K2.O2*L*w*h
                                                   *(C.O2.sat-C.O2))), # gas exchange
       processes        = list(gro.ALG.NH4,gro.ALG.NO3,resp.ALG, #   
                               nitri1,nitri2,miner.POM))

reach3 <- 
   new(Class            = "reactor",
       name             = "R3",
       volume.ini       = expression(L*w*h),
       area             = expression(L*w),
       conc.pervol.ini  = list(C.HPO4 = expression(C.HPO4.ini),  # gP/m3
                               C.NH4  = expression(C.NH4.ini),   # gN/m3
                               C.NO2  = expression(C.NO2.ini),   # gN/m3
                               C.NO3  = expression(C.NO3.ini),   # gN/m3
                               C.O2   = expression(C.O2.ini)),   # gN/m3
       input            = list(C.O2   = expression(K2.O2*L*w*h
                                                   *(C.O2.sat-C.O2))), # gas ex.
       outflow          = expression(Q.in*86400),
       processes        = list(gro.ALG.NH4,gro.ALG.NO3,resp.ALG, #  
                               nitri1,nitri2,miner.POM))


# Definition of links:
# ====================

reach1.reach2 <-
   new(Class      = "link",
       name       = "R1R2",
       from       = "R1",
       to         = "R2",
       flow       = expression(Q.in*86400))

reach2.reach3 <-
   new(Class      = "link",
       name       = "R2R3",
       from       = "R2",
       to         = "R3",
       flow       = expression(Q.in*86400))


# Definition of system:
# =====================

# River system:
# -------------

system <- new(Class    = "system",
              name     = "River",
              reactors = list(reach1,reach2,reach3),
              links    = list(reach1.reach2,reach2.reach3),
              cond     = cond,
              param    = param,
              t.out    = seq(0,3,by=0.02))


# Perform simulation:
# ===================

res <- calcres(system)


# Plot results:
# =============
                 
#plotres(res)              # plot to screen

#plotres(res,colnames=list(c("C.NH4.R1", "C.NH4.R2", "C.NH4.R3"),
#                          c("C.NO2.R1", "C.NO2.R2", "C.NO2.R3"),
#                          c("C.NO3.R1", "C.NO3.R2", "C.NO3.R3"),
#                          c("C.HPO4.R1","C.HPO4.R2","C.HPO4.R3"),
#                          c("C.O2.R1",  "C.O2.R2",  "C.O2.R3")))

plotres(res,colnames=list(c("C.NH4.R1", "C.NH4.R2", "C.NH4.R3"),
                          c("C.NO2.R1", "C.NO2.R2", "C.NO2.R3"),
                          c("C.NO3.R1", "C.NO3.R2", "C.NO3.R3"),
                          c("C.HPO4.R1","C.HPO4.R2","C.HPO4.R3"),
                          c("C.O2.R1",  "C.O2.R2",  "C.O2.R3")),
         file     = "rivermodel_simple.pdf",
         width    = 8,
         height   = 6)

