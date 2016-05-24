# ==============================================================================
# Simple Lake Model
# ==============================================================================


# Load package:
# =============

library(ecosim)


# Definition of parameters:
# =========================


param    <- list(k.gro.ALG   = 1,        # 1/d
                 k.gro.ZOO   = 0.8,      # m3/gDM/d
                 k.death.ALG = 0.4,      # 1/d
                 k.death.ZOO = 0.08,     # 1/d
                 K.HPO4      = 0.002,    # gP/m3
                 Y.ZOO       = 0.2,      # gDM/gDM
                 alpha.P.ALG = 0.002,    # gP/gDM
                 A           = 8.5e+006, # m2
                 h.epi       = 4,        # m
                 Q.in        = 4,        # m3/s
                 C.ALG.ini   = 0.05,     # gDM/m3
                 C.ZOO.ini   = 0.1,      # gDM/m3
                 C.HPO4.ini  = 0.02,     # gP/m3
                 C.HPO4.in   = 0.04)     # gP/m3             


# Definition of transformation processes:
# =======================================

# Growth of algae:
# ----------------

gro.ALG   <- new(Class  = "process",
                 name   = "Growth of algae",
                 rate   = expression(k.gro.ALG
                                     *C.HPO4/(K.HPO4+C.HPO4)
                                     *C.ALG),
                 stoich = list(C.ALG  = expression(1),              # gDM/gDM
                               C.HPO4 = expression(-alpha.P.ALG)))  # gP/gDM

# Death of algae:
# ---------------

death.ALG <- new(Class = "process",
                 name   = "Death of algae",
                 rate   = expression(k.death.ALG*C.ALG),
                 stoich = list(C.ALG  = expression(-1)))            # gDM/gDM

# Growth of zooplankton:
# ----------------------

gro.ZOO   <- new(Class  = "process",
                 name   = "Growth of zooplankton",
                 rate   = expression(k.gro.ZOO
                                     *C.ALG
                                     *C.ZOO),
                 stoich = list(C.ZOO  = expression(1),              # gDM/gDM
                               C.ALG  = expression(-1/Y.ZOO)))      # gP/gDM

# Death of zooplankton:
# ---------------------

death.ZOO <- new(Class  = "process",
                 name   = "Death of zooplankton",
                 rate   = expression(k.death.ZOO*C.ZOO),
                 stoich = list(C.ZOO  = expression(-1)))            # gDM/gDM


# Definition of reactor:
# ======================

# Epilimnion:
# -----------

epilimnion <- 
   new(Class            = "reactor",
       name             = "Epilimnion",
       volume.ini       = expression(A*h.epi),
       conc.pervol.ini  = list(C.HPO4 = expression(C.HPO4.ini),     # gP/m3
                               C.ALG  = expression(C.ALG.ini),      # gDM/m3
                               C.ZOO  = expression(C.ZOO.ini)),     # gDM/m3
       inflow           = expression(Q.in*86400),                   # m3/d
       inflow.conc      = list(C.HPO4 = expression(C.HPO4.in),
                               C.ALG  = 0,
                               C.ZOO  = 0),
       outflow          = expression(Q.in*86400),
       processes        = list(gro.ALG,death.ALG,gro.ZOO,death.ZOO))


# Definition of system:
# =====================

# Lake system:
# ------------

system <- new(Class    = "system",
              name     = "Lake",
              reactors = list(epilimnion),
              param      = param,
              t.out    = seq(0,365,by=1))


# Perform simulation:
# ===================

res1 <- calcres(system)


# Plot results:
# =============
                 
#plotres(res1)              # plot to screen

plotres(res=res1,file="lakemodel_simple_11.pdf")  # plot to pdf file

plotres(res=res1, colnames=c("C.ALG", "C.ZOO"))

#plotres(res=res1, colnames=list("C.HPO4",c("C.ALG", "C.ZOO")))

#plotres(res=res1[1:100,], colnames=list("C.HPO4",c("C.ALG", "C.ZOO")))

plotres(res      = res1,    # plot to pdf file
        colnames = list("C.HPO4",c("C.ALG","C.ZOO")),
        file     = "lakemodel_simple_12.pdf",
        width    = 8,
        height   = 4)



# ==============================================================================
# Model Extension
# ==============================================================================


# Growth of algae extended process:
# ---------------------------------

gro.ALG.ext <-
   new(Class  = "process",
       name   = "Growth of algae extended",
       rate   = expression(k.gro.ALG
                           *exp(beta.ALG*(T-T0))
                           *C.HPO4/(K.HPO4+C.HPO4)
                           *log((K.I+I0)
                                /(K.I+I0*exp(-(lambda.1+lambda.2*C.ALG)*h.epi)))
                            /((lambda.1+lambda.2*C.ALG)*h.epi)
                           *C.ALG),
       stoich = list(C.ALG  = 1,                            # gDM/gDM
                     C.HPO4 = expression(-alpha.P.ALG)))    # gP/gDM

# Change processes in the reactor "epilimnion"
epilimnion@processes <- list(gro.ALG.ext,death.ALG,gro.ZOO,death.ZOO)


# Changed parameter values and additional parameters:
# ===================================================

param    <- list(k.gro.ALG   = 1.4,      # 1/d
                 k.gro.ZOO   = 0.6,      # m3/gDM/d
                 k.death.ALG = 0.4,      # 1/d
                 k.death.ZOO = 0.08,     # 1/d
                 K.HPO4      = 0.002,    # gP/m3
                 Y.ZOO       = 0.2,      # gDM/gDM
                 alpha.P.ALG = 0.002,    # gP/gDM
                 A           = 8.5e+006, # m2
                 h.epi       = 4,        # m
                 Q.in        = 4,        # m3/s
                 C.ALG.ini   = 0.05,     # gDM/m3
                 C.ZOO.ini   = 0.1,      # gDM/m3
                 C.HPO4.ini  = 0.02,     # gP/m3
                 C.HPO4.in   = 0.04,     # gP/m3
                 beta.ALG    = 0.046,    # 1/degC
                 T0          = 20,       # degC
                 K.I         = 30,       # W/m2
                 lambda.1    = 0.05,     # 1/m
                 lambda.2    = 0.03)     # m2/gDM                  


# Add additional parameters for environmental conditions:
# =======================================================

param <- c(param,
           list(t.max  = 230,
                I0.min = 25,
                I0.max = 225,
                T.min  = 5,
                T.max  = 25))
  
# Change the parameters in the system         
system@param <- param

# Change the environmental conditions in the epilimnion

epilimnion@cond <- list(I0 = expression(0.5*(I0.min+I0.max)+
                                        0.5*(I0.max-I0.min)*
                                        cos(2*pi/365.25*(t-t.max))),   # W/m2
                        T  = expression(0.5*(T.min+T.max)+
                                        0.5*(T.max-T.min)*
                                         cos(2*pi/365.25*(t-t.max))))  # degC
system@reactors <- list(epilimnion)

t.out    = seq(0,1461,by=1)
system@t.out <- t.out

res2 <- calcres(system)
       
                      
#plotres(res2)                       # plot to screen

plotres(res=res2,file="lakemodel_simple_21.pdf")  # plot to pdf file    

plotres(res=res2, colnames=c("C.ALG", "C.ZOO"))

#plotres(res=res2, colnames=list("C.HPO4",c("C.ALG", "C.ZOO")))

#plotres(res=res2[1:365,], colnames=list("C.HPO4",c("C.ALG", "C.ZOO")))

plotres(res      = res2,    # plot to pdf file
        colnames = list("C.HPO4",c("C.ALG","C.ZOO")),
        file     = "lakemodel_simple_22.pdf",
        width    = 8,
        height   = 4)

