################################################################################
#  pH model of the Scheldt estuary                                             #
#  Hofmann AF, Meysman FJR, Soetaert K, Middelburg J, 2008.                    #
#  A step-by-step procedure for pH model construction in aquatic systems       #
#  Biogeosciences 5, 227-251                                                   #
#                                                                              #
#                               STEP 2 - FNA                                   #
#  Full numerical approach - pH model written as a set of                      #
#  differential algebraic equations, solved with DAE solver daspk              #
#  Implementation: Andreas Hofmann, Karline Soetaert - NIOZ                    #
################################################################################

# load parameters, dissociation constants, initial conditions,
# the model transport function and function TA_estimate, to estimate alkalinity
# Do make sure that this file is in the working directory (or use setwd(""))


source('Schelde_pars.R')


################################################################################
#                   DIFFERENTIAL ALGEBRAIC EQUATIONS                           #
################################################################################

FNAResidual <- function (tt, state, dy, parms, scenario = "B1") {
  with (as.list(c(state, dy, parms)), {
    pH      <- -log10(H*1e-6)
    TA      <- HCO3 + 2*CO3 + NH3 - H
    SumCO2  <- CO2 + HCO3 + CO3
    SumNH4  <- NH4 + NH3

    #--------------------------
    # PHYSICAL PROCESSES
    #--------------------------

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
           AddNH4NO3 <- SpillNH4NO3   # NH4+NO3- - tanker addition
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
    # RESIDUALS OF RATE OF CHANGES
    #--------------------------
    # 9 unknowns (dOM,dO2,dNO3,dCO2,dHCO3,dCO3,dNH4,dNH3,dH) - 9 equations
    # of simple state variables
    ROM  <- - dOM  + TOM         - ROx
    RO2  <- - dO2  + TO2  + EO2  - ROxCarbon - 2*RNit
    RNO3 <- - dNO3 + TNO3                    +   RNit  + AddNH4NO3

    # of summed quantities
    RSumCO2 <- -dCO2 -dHCO3 -dCO3 + TSumCO2 + ECO2 + ROxCarbon
    RSumNH4 <- -dNH3 -dNH4        + TSumNH4 + ENH3 + ROx - RNit + AddNH3 + AddNH4NO3
    RTA     <- -dHCO3-2*dCO3-dNH3 +dH + TTA + ENH3 + ROx - 2*RNit + AddNH3

    # algebraic equations: equilibrium equations
    EquiCO2 <- H * HCO3 - K1CO2 * CO2
    EquiHCO3<- H * CO3 - K2CO2  * HCO3
    EquiNH4 <- H * NH3 - KNH4   * NH4

    #--------------------------
    # Output variables: The pH, alkalinity and other summed quantities
    #--------------------------
    return(list(c(ROM, RO2, RNO3, RSumCO2, RSumNH4, RTA, EquiCO2, EquiHCO3, EquiNH4),
    c(pH = pH, TA = TA, SumCO2 = SumCO2, SumNH4 = SumNH4)))
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
# initial conditions
#---------------------

   H_ini <- 10^-pH_ini * 1e6    
       H <- H_ini
 NH3_ini <- KNH4/(KNH4+H)*SumNH4_ini
 NH4_ini <- SumNH4_ini - NH3_ini
 CO2_ini <- H*H/(H*K1CO2 + H*H + K1CO2*K2CO2)*SumCO2_ini
HCO3_ini <- H*K1CO2/(H*K1CO2 + H*H + K1CO2*K2CO2)*SumCO2_ini
 CO3_ini <- K1CO2*K2CO2/(H*K1CO2 + H*H + K1CO2*K2CO2)*SumCO2_ini
  TA_ini <- HCO3_ini + 2*CO3_ini + NH3_ini - H_ini

# Initial conditions for the state variables AND their rates of change
y <- c(OM = OM_ini, O2 = O2_ini, NO3 = NO3_ini, H = H_ini,
       NH4 = NH4_ini, NH3 = NH3_ini, CO2 = CO2_ini, HCO3 = HCO3_ini, CO3 = CO3_ini)
dy <- c(dOM = 0, dO2 = 0, dNO3 = 0, dH = 0, dNH4 = 0, dNH3 = 0, dCO2 = 0, dHCO3 = 0, dCO3 = 0)

#---------------------
# run the model
#---------------------

times <- c(0, 350:405)

outA <- daspk(y = y, times, res = FNAResidual, dy = dy, nalg = 3, 
  parms = phPars, scenario = "A", hmax = 1)
outB <- daspk(y = y, times, res = FNAResidual, dy = dy, nalg = 3, 
  parms = phPars, scenario = "B1", hmax = 1)
outC <- daspk(y = y, times, res = FNAResidual, dy = dy, nalg = 3, 
  parms = phPars, scenario = "C",  hmax = 1)

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
