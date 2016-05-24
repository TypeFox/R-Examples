################################################################################
#  pH model of the Scheldt estuary                                             #
#  Hofmann AF, Meysman FJR, Soetaert K, Middelburg J, 2008.                    #
#  A step-by-step procedure for pH model construction in aquatic systems       #
#  Biogeosciences                                                              #
#                                                                              #
#  MODEL PARAMETERS, INITIAL CONDITIONS, COMMON MODEL ROUTINES                 #
#  Implementation: Andreas Hofmann, Karline Soetaert - NIOZ                    #
################################################################################

require(deSolve)

################################################################################
# Global Physical parameters ##
################################################################################
Q	     <- 8640000      # m3/d	discharge
V      <- 108798000    # m3	volume
Eprime <- 13824000     # m3/d	averaged bulk-dispersion coefficient, 160 m3/s)

################################################################################
# boundary conditions
################################################################################
# upper boundary
OM_up	    <- 50		     # umol/kg-soln
NO3_up    <- 350		   # umol/kg-soln
O2_up	    <- 70		     # umol/kg-soln
pH_up     <- 7.6
SumNH4_up <- 80		     # umol/kg-soln
SumCO2_up <- 7100      # umol/kg-soln


# lower boundary - pH and alkalinity are consistent
OM_down	    <- 25		   # umol/kg-soln
NO3_down    <- 260		 # umol/kg-soln
O2_down	    <- 240		 # umol/kg-soln
pH_down     <- 7.92
SumNH4_down <- 7		   # umol/kg-soln
SumCO2_down <- 4400	   # umol/kg-soln

################################################################################
# initial conditions: as derived from steady state run; pH and alkinity consistent
################################################################################

OM_ini      <- 31.9688	  # umol/kg-soln
NO3_ini     <- 340.235	  # umol/kg-soln
O2_ini      <- 157.922	  # umol/kg-soln
pH_ini      <- 7.7        #
SumNH4_ini  <- 35.8406	  # umol/kg-soln
SumCO2_ini  <- 6017.28	  # umol/kg-soln

################################################################################
#                               MODEL PARAMETERS                               #
################################################################################

phPars <- c(
  KL	      = 0.28  ,  #   1/d	            proportionality factor for air-water exchange
  rOM	      = 0.1   ,  #   1/d	            first-order oxic mineralisation rate of organic matter
  rNitri    = 0.26  ,  #   1/d              first order nitrification rate (with resp. to Ammonium)
  ksO2      = 20.0  ,  #   umol-O2/kg-soln  monod half-saturation constant Oxygen (ox min & nit)
  k1        = 1e3   ,  #   1/d              "instantaneous" rate for forward equilibrium reactions
  C_Nratio  = 8     ,  #   mol C/mol N      C:N ratio oforganic matter
  rDenit	  = 0.2   ,  #   1/d		 first    order mineralisation due to denit rate (w.r.t. OM)
  ksNO3	    = 45	   , #   umol-NO3/kg	    monod half-saturation constant nitrate denitrification
  ksO2inhib = 22	   , #   umol-02/kg	      monod inhibition term oxygen

  # saturated concentrations - calculated for T=12 and S=5 #
  CO2sat	  = 19     , #   umol/kg-soln
  O2sat	    = 325    , #   umol/kg-soln
  NH3sat    = 0.0001 , #   umol/kg-soln

  ################################################################################
  ## DIFFERENT SCENARIOS:
  #	@ A	 decreased waste load due to a sewage treatement plant in Brussels
  #	@ B1	 a 10000 ton fertilizer (NH4+/NO3-) ship sinks: different modelling approach (extra NH4NO3 addition)
  #	@ C	 a 10000 ton NH3 ship sinks: modelling approach 1 (extra NH3 addition)
  ################################################################################

  # Scenario A: Brussels wastewater treatment plant scenario reduces upstream conc of OM #
  OM_up_A	    = 25	 , # umol/kg-soln

  # Scenario B1: Ammonium-Nitrate (fertilizer) tank ship scenario: #
  # model it as extra NH4+ and NO3 - addition of 10000 tpms#
  SpillNH4NO3 = ((10000 * 1000000)/(18 + 62)) * # Total substance in mol over 10 days
                 1000000 / (V * 1000) / 10,              # Conc in umol/kg per day

  # Scenario C: NH3 (Ammonia) tank ship scenario (10000 tons NH3 input)  #
  SpillNH3    = ((10000 * 1000000) / 17) * # Total substance in mol/10 days
                   1000000 / (V * 1000) / 10     # Conc in umol/kg per day
)


################################################################################
# Dissociation constants
################################################################################

require(seacarb)
# Temperature, salinity settings
Temp      <- 12      # dg C
Sal       <- 5       #

K1CO2     <- K1(S = Sal, T = Temp, P = 0)*1e6  # umol/kg-soln
K2CO2     <- K2(S = Sal, T = Temp, P = 0)*1e6  # umol/kg-soln
KNH4      <- Kn(S = Sal, T = Temp, P = 0)*1e6  # umol/kg-soln
KW        <- Kw(S = Sal, T = Temp, P = 0)*1e12 # (mol/kg-soln)^2


################################################################################
#                           COMMON MODEL FUNCTIONS                             #
################################################################################

# Advective-dispersive transport
Transport <- function (y, y.up, y.down) {  # Q: discharge, m3/d; Eprime: bulk dispersion coefficient, V: Volume
  Input <- Q * c(y.up, y) - Eprime * diff(c(y.up, y, y.down))
  dy    <- -diff(Input)/V
  return(dy)
}

# Estimate alkalinity based on pH, sum CO2, sum NH4 
TA_estimate <- function(pH, DIC, SumNH4) {
  H <- 10^(-pH)*1e6
  HCO3 <- H*K1CO2/(H*K1CO2 + H*H + K1CO2*K2CO2)*DIC
  CO3 <- K1CO2*K2CO2/(H*K1CO2 + H*H + K1CO2*K2CO2)*DIC
  NH3 <- KNH4/(KNH4+H)*SumNH4
  return(as.double(HCO3 + 2*CO3 + NH3 - H))   # Total alkalinity
}

