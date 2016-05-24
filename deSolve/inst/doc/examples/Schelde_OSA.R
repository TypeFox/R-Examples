################################################################################
#  pH model of the Scheldt estuary                                             #
#  Hofmann AF, Meysman FJR, Soetaert K, Middelburg J, 2008.                    #
#  A step-by-step procedure for pH model construction in aquatic systems       #
#  Biogeosciences 5, 227-251                                                   #
#                                                                              #
#                               STEP 3 -OSA                                    #
#  Operator splitter approach - pH model written as a set of                   #
#  ordinary differential equations, solved with ODE solver vode                #
#  Each time step the pH is solved at equilibrium, using uniroot               #   
#  Implementation: Andreas Hofmann, Karline Soetaert - NIOZ                    #     
################################################################################

# load parameters, dissociation constants, initial conditions,
# the model transport function and function TA_estimate, to estimate alkalinity
# Do make sure that this file is in the working directory (or use setwd(""))

source('Schelde_pars.R')

################################################################################
#                                  UTILITIES                                   #
################################################################################

# Function that estimates discrepancy between estimated and true total alkalinity 
# Root of this function = solution of equilibrium pH

pHfunction <- function(pH, DIC, TA, SumNH4)  return(TA-TA_estimate(pH, DIC, SumNH4))

################################################################################
#                     ORDINARY DIFFERENTIAL EQUATIONS                          #
################################################################################

OSAmodel <- function (tt, state, parms, scenario="B1") {
  with (as.list(c(state, parms)), {

    pH <- uniroot (pHfunction, interval = c(6, 10), tol=1e-20,
                   DIC=SumCO2, TA=TA, SumNH4=SumNH4)$root

    #--------------------------
    # PHYSICAL PROCESSES
    #--------------------------
    H   <- 10^(-pH) * 1e6
    CO2 <- H*H/(H*K1CO2 + H*H + K1CO2*K2CO2)*SumCO2

    NH3 <- KNH4/(KNH4+H)*SumNH4
    NH4 <- SumNH4 - NH3

    # air-water exchange
    ECO2    <- KL * (CO2sat - CO2)
    EO2     <- KL * (O2sat  - O2)
    ENH3    <- KL * (NH3sat - NH3)

    # Transport
    TO2     <- Transport(O2,     O2_up,     O2_down)
    TNO3    <- Transport(NO3,    NO3_up,    NO3_down)
    TTA     <- Transport(TA,     TA_up,     TA_down)
    TSumCO2 <- Transport(SumCO2, SumCO2_up, SumCO2_down)
    TSumNH4 <- Transport(SumNH4, SumNH4_up, SumNH4_down)

    # Wastewater treatment plant in Brussels scenario
    if (scenario == "A" && tt > 365) {
           TOM <- Transport(OM,     OM_up_A, OM_down)
    } else TOM <- Transport(OM,     OM_up  , OM_down)

    # Spills
    if (scenario == "B1" &&  (tt > 360 && tt < 370)) {
           AddNH4NO3 <- SpillNH4NO3     # NH4+NO3- - tanker addition
    } else AddNH4NO3 <- 0

    if (scenario == "C"  && (tt > 360 && tt < 370))  {
           AddNH3    <- SpillNH3        # NH3 - tanker input
    } else AddNH3    <- 0


    #--------------------------
    # BIOGEOCHEMICAL PROCESSES:
    #--------------------------

    # Oxic mineralisation
    ROx       <- rOM * OM * (O2/(O2 + ksO2))
    ROxCarbon <- ROx * C_Nratio

    # Nitrification
    RNit  <- rNitri * NH4 * (O2/(O2 + ksO2))

    #--------------------------
    # RATE OF CHANGE
    #--------------------------

    dOM     <-  TOM         - ROx
    dO2     <-  TO2  + EO2  - ROxCarbon - 2*RNit
    dNO3    <-  TNO3                    +   RNit  + AddNH4NO3
    dSumCO2 <-  TSumCO2 + ECO2 + ROxCarbon
    dSumNH4 <-  TSumNH4 + ENH3 + ROx - RNit + AddNH3 + AddNH4NO3
    dTA     <-  TTA  + ENH3    + ROx-2*RNit + AddNH3

    return(list(c(dOM, dO2, dNO3, dTA, dSumNH4, dSumCO2),
    c(pH=pH, CO2=CO2, NH3=NH3, NH4=SumNH4-NH3)))
  })
}




################################################################################
#                             MODEL APPLICATIONS                               #
################################################################################

#---------------------
# Akalinity at boundaries
#---------------------

TA_down<- TA_estimate(pH_down, SumCO2_down, SumNH4_down)
TA_up  <- TA_estimate(pH_up  , SumCO2_up  , SumNH4_up)


#---------------------
# initial conditions
#---------------------
TA_ini <- TA_estimate(pH_ini , SumCO2_ini , SumNH4_ini)

state <- c(OM=OM_ini, O2=O2_ini, NO3=NO3_ini, TA=TA_ini,
           SumNH4=SumNH4_ini, SumCO2=SumCO2_ini)

#---------------------
# run model
#---------------------
times <- c(0, 350:405)

outA <- vode(state, times, OSAmodel, phPars, scenario = "A" ,  hmax = 1)
outB <- vode(state, times, OSAmodel, phPars, scenario = "B1",  hmax = 1)
outC <- vode(state, times, OSAmodel, phPars, scenario = "C" ,  hmax = 1)

#---------------------
# plot model output
#---------------------
par(mfrow = c(3, 4), mar = c(1, 2, 2, 1),  oma = c(3, 3, 3, 0))
Selection <- c("pH","TA","SumCO2","O2")
plot(outA, mfrow = NULL, xlim = c(350,405), type = "l",
  xaxt = "n", which = Selection)
plot(outB, mfrow = NULL, xlim = c(350,405), type = "l",
  xaxt = "n", which = Selection)
plot(outC, mfrow = NULL, xlim = c(350,405), type = "l",
  xaxt = "n", which = Selection)

mtext(side = 1, outer = TRUE, "time, d", line = 2, cex = 1.2)
mtext(side = 2, at = 0.2, outer = TRUE, "Scenario C", line = 1.5, cex = 1.2)
mtext(side = 2, at = 0.5, outer = TRUE, "Scenario B", line = 1.5, cex = 1.2)
mtext(side = 2, at = 0.8, outer = TRUE, "Scenario A", line = 1.5, cex = 1.2)

mtext(side = 3, at = 0.125, outer = TRUE, "pH, -", line = 1, cex = 1.2)
mtext(side = 3, at = 0.375, outer = TRUE, "TA, µmol/kg", line = 1, cex = 1.2)
mtext(side = 3, at = 1-0.375, outer = TRUE, "CO2, µmol/kg", line = 1, cex = 1.2)
mtext(side = 3, at = 1-0.125, outer = TRUE, "O2, µmol/kg", line = 1, cex = 1.2)
