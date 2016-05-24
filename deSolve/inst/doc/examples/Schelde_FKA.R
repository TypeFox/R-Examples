################################################################################
#  pH model of the Scheldt estuary                                             #
#  Hofmann AF, Meysman FJR, Soetaert K, Middelburg J, 2008.                    #
#  A step-by-step procedure for pH model construction in aquatic systems       #
#  Biogeosciences 5, 227-251                                                   #
#                                                                              #
#                               STEP 1 - FKA                                   #
#  Full kinetic approach - pH model written as a set of stiff                  #
#  ordinary differential equations, solved with ODE solver vode                #
#  Implementation: Andreas Hofmann, Karline Soetaert - NIOZ                    #
################################################################################

# load parameters, dissociation constants, initial conditions,
# the model transport function and function TA_estimate, to estimate alkalinity
# Do make sure that this file is in the working directory (or use setwd(""))

source('Schelde_pars.R')

################################################################################
#                               MODEL EQUATIONS                                #
################################################################################

FKAmodel <- function (tt, state, parms, scenario="B1") {
  with (as.list(c(state, parms)), {

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

    TH      <- Transport(H,      H_up,      H_down)
    TCO2    <- Transport(CO2,    CO2_up,    CO2_down)
    THCO3   <- Transport(HCO3,   HCO3_up,   HCO3_down)
    TCO3    <- Transport(CO3,    CO3_up,    CO3_down)

    TNH3    <- Transport(NH3,    NH3_up,    NH3_down)
    TNH4    <- Transport(NH4,    NH4_up,    NH4_down)

    # Wastewater treatment plant in Brussels scenario
    if (scenario == "A" && tt > 365) {
           TOM <- Transport(OM,     OM_up_A, OM_down)
    } else TOM <- Transport(OM,     OM_up  , OM_down)

    # Spills
    if (scenario == "B1" &&  (tt > 360 && tt < 370)) {
           AddNH4NO3 <- SpillNH4NO3    # NH4+NO3- - tanker addition
    } else AddNH4NO3 <- 0

    if (scenario == "C"  && (tt > 360 && tt < 370))  {
           AddNH3    <- SpillNH3       # NH3 - tanker input
    } else AddNH3    <- 0


    #--------------------------
    # BIOGEOCHEMICAL PROCESSES:
    #--------------------------

    # Oxic mineralisation
    ROx       <- rOM * OM * (O2/(O2 + ksO2))
    ROxCarbon <- ROx * C_Nratio

    # Nitrification
    RNit  <- rNitri * NH4 * (O2/(O2 + ksO2))

    # "equilibrium reactions": k1 arbitrarily high
    RCO2  <- k1*CO2  - k1/K1CO2* H * HCO3
    RHCO3 <- k1*HCO3 - k1/K2CO2* H * CO3
    RNH4  <- k1*NH4  - k1/KNH4 * H * NH3

    #--------------------------
    # RATE OF CHANGE
    #--------------------------

    dOM     <-  TOM         - ROx
    dO2     <-  TO2  + EO2  - ROxCarbon - 2*RNit
    dNO3    <-  TNO3                    +   RNit  + AddNH4NO3

    dCO2    <-  TCO2 + ECO2 + ROxCarbon - RCO2
    dHCO3   <-  THCO3                   + RCO2 - RHCO3
    dCO3    <-  TCO3                    + RHCO3

    dNH3    <-  TNH3 + ENH3 + ROx       + RNH4 + AddNH3
    dNH4    <-  TNH4 - RNit             - RNH4 + AddNH4NO3

    dH      <-  TH   + 2*RNit + RCO2    + RHCO3 + RNH4

    #--------------------------
    # Output variables: The pH, alkalinity and other summed quantities
    #--------------------------
    pH      <- -log10(H*1e-6)
    TA      <- HCO3 + 2*CO3 + NH3 - H
    SumCO2  <- CO2 + HCO3 + CO3
    SumNH4  <- NH4 + NH3
    return(list(c(dOM, dO2, dNO3, dH, dNH4, dNH3, dCO2, dHCO3, dCO3),
    c(pH=pH, TA=TA, SumCO2=SumCO2, SumNH4=SumNH4)))
  })
}




################################################################################
#                             MODEL APPLICATIONS                               #
################################################################################

#---------------------
# Extra Boundary conditions
#---------------------

# The speciation of DIC and sum(ammonium),  calculated consistently with pH_up
    H_up <- 10^-pH_up * 1e6    # umol/kg solution
       H <- H_up
  NH3_up <- KNH4/(KNH4+H)*SumNH4_up
  NH4_up <- SumNH4_up - NH3_up
  CO2_up <- H*H/(H*K1CO2 + H*H + K1CO2*K2CO2)*SumCO2_up
 HCO3_up <- H*K1CO2/(H*K1CO2 + H*H + K1CO2*K2CO2)*SumCO2_up
  CO3_up <- K1CO2*K2CO2/(H*K1CO2 + H*H + K1CO2*K2CO2)*SumCO2_up


# calculated consistently with pH_down:
   H_down <- 10^-pH_down * 1e6    # umol/kg solution
        H <- H_down
 NH3_down <- KNH4/(KNH4+H)*SumNH4_down
 NH4_down <- SumNH4_down - NH3_down
 CO2_down <- H*H/(H*K1CO2 + H*H + K1CO2*K2CO2)*SumCO2_down
HCO3_down <- H*K1CO2/(H*K1CO2 + H*H + K1CO2*K2CO2)*SumCO2_down
 CO3_down <- K1CO2*K2CO2/(H*K1CO2 + H*H + K1CO2*K2CO2)*SumCO2_down

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

state <- c(OM=OM_ini, O2=O2_ini, NO3=NO3_ini, H=H_ini,
           NH4=NH4_ini, NH3=NH3_ini, CO2=CO2_ini, HCO3=HCO3_ini, CO3=CO3_ini)

#---------------------
# run model
#---------------------

times <- c(0, 350:405)

outA <- vode(state, times, FKAmodel, phPars, scenario = "A" , hmax = 1)
outB <- vode(state, times, FKAmodel, phPars, scenario = "B1", hmax = 1)
outC <- vode(state, times, FKAmodel, phPars, scenario = "C" , hmax = 1)

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
