################################################################################
#  pH model of the Scheldt estuary                                             #
#  Hofmann AF, Meysman FJR, Soetaert K, Middelburg J, 2008.                    #
#  A step-by-step procedure for pH model construction in aquatic systems       #
#  Biogeosciences 5, 227-251                                                   #
#                                                                              #
#                               STEP 4 -DSA                                    #
#  Direct substitution approach - pH model written as a set of                 #
#  ordinary differential equations, solved with ODE solver vode                #
#  Hplus is a state variable; the model is not stiff                           #        
#  Implementation: Andreas Hofmann, Karline Soetaert - NIOZ                    #
################################################################################

# load parameters, dissociation constants, initial conditions,
# the model transport function and function TA_estimate, to estimate alkalinity
# Do make sure that this file is in the working directory (or use setwd(""))

source('Schelde_pars.R')


################################################################################
#                     ORDINARY DIFFERENTIAL EQUATIONS                          #
################################################################################

DSAmodel <- function (tt, state, parms, scenario = "B1") {
  with (as.list(c(state, parms)), {

    pH <- -log10(H*1e-6)
    TA <- TA_estimate(pH, SumCO2, SumNH4)

    #--------------------------
    # PHYSICAL PROCESSES
    #--------------------------
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
           AddNH4NO3 <- SpillNH4NO3    # NH4+NO3- - tanker addition
    } else AddNH4NO3 <- 0

    if (scenario == "C"  && (tt > 360 && tt < 370))  {
           AddNH3    <- SpillNH3      # NH3 - tanker input:
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

    # rate of change of pH:
    dTAdSumCO2  <- (H*K1CO2 + (2*K1CO2*K2CO2))/((H*K1CO2) + (K1CO2*K2CO2) + (H*H))
    dTAdSumNH4  <- KNH4 / (KNH4 + H)

    dHCO3dH  <- ((K1CO2/((H*K1CO2) + (K1CO2*K2CO2) +  (H*H))) -
       ((H*K1CO2*((2*H)+K1CO2))/(((H*K1CO2) + (K1CO2*K2CO2) +  (H*H))^2)))* SumCO2
    dCO3dH   <- -((K1CO2*K2CO2*((2*H)+K1CO2))/
                 (((H*K1CO2) + (K1CO2*K2CO2) +  (H*H))^2)) * SumCO2
    dNH3dH   <- -(KNH4  / ((H*H)+(2*H*KNH4)+(KNH4*KNH4)))  * SumNH4

    dHdH     <- 1
    dTAdH    <- dHCO3dH + 2*dCO3dH + dNH3dH - dHdH
    dH       <- ((ROx - 2*RNit + ENH3 + AddNH3 + TTA) -
         ((dTAdSumCO2*dSumCO2) + (dTAdSumNH4*dSumNH4)))/dTAdH

    return(list(c(dOM, dO2, dNO3, dH, dSumNH4, dSumCO2),
               c(TA=TA, pH=pH, CO2=CO2, NH3=NH3, NH4=SumNH4-NH3)))
  })
}


################################################################################
#                             MODEL APPLICATIONS                               #
################################################################################

#---------------------
# alkalinity at the boundaries
#---------------------

TA_up    <- TA_estimate(pH_up, SumCO2_up, SumNH4_up)
TA_down  <- TA_estimate(pH_down, SumCO2_down, SumNH4_down)

#---------------------
# the initial conditions
#---------------------
H_ini <- 10^(-pH_ini)*1e6

state <- c(OM=OM_ini, O2=O2_ini, NO3=NO3_ini, H=H_ini,
           SumNH4=SumNH4_ini, SumCO2=SumCO2_ini)

#---------------------
# run model - three scenarios
#---------------------
times <- c(0, 350:405)

outA <- vode(state, times, DSAmodel, phPars, scenario = "A",  hmax = 1)
outB <- vode(state, times, DSAmodel, phPars, scenario = "B1", hmax = 1)
outC <- vode(state, times, DSAmodel, phPars, scenario = "C" , hmax = 1)

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
