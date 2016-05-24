.BROScalc <- function(FUELTYPE, FFMC, BUI, WSV, FMC, SFC, PC, PDF, CC, CBH){
  #############################################################################
  # Description:
  #   Calculate the Back Fire Spread Rate. 
  #
  #   All variables names are laid out in the same manner as Forestry Canada 
  #   Fire Danger Group (FCFDG) (1992). Development and Structure of the 
  #   Canadian Forest Fire Behavior Prediction System." Technical Report 
  #   ST-X-3, Forestry Canada, Ottawa, Ontario.
  #
  # Args:
  #   FUELTYPE: The Fire Behaviour Prediction FuelType
  #   FFMC:     Fine Fuel Moisture Code
  #   BUI:      Buildup Index
  #   WSV:      Wind Speed Vector
  #   FMC:      Foliar Moisture Content
  #   SFC:      Surface Fuel Consumption
  #   PC:       Percent Conifer
  #   PDF:      Percent Dead Balsam Fir
  #   CC:       Degree of Curing (just "C" in FCFDG 1992)
  #   CBH:      Crown Base Height
  
  # Returns:
  #   BROS:     Back Fire Spread Rate
  #
  #############################################################################
  #Eq. 46 (FCFDG 1992)
  #Calculate the FFMC function from the ISI equation
  m <- 147.2 * (101 - FFMC) / (59.5 + FFMC)
  #Eq. 45 (FCFDG 1992)
  fF <- 91.9 * exp(-0.1386 * m) * (1.0 + (m**5.31) / 4.93e7)
  #Eq. 75 (FCFDG 1992)
  #Calculate the Back fire wind function
  BfW <- exp(-0.05039 * WSV)
  #Calculate the ISI associated with the back fire spread rate
  #Eq. 76 (FCFDG 1992)
  BISI <- 0.208 * BfW * fF
  #Eq. 77 (FCFDG 1992)
  #Calculate final Back fire spread rate
  BROS <- .ROScalc(FUELTYPE, BISI, BUI, FMC, SFC, PC, PDF, CC, CBH)
  
  return(BROS)
}