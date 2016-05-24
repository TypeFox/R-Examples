.DISTtcalc <- function(FUELTYPE, ROSeq, HR, CFB) {
  #############################################################################
  # Description:
  #   Calculate the Head fire spread distance at time t. In the documentation
  #   this variable is just "D".
  #
  #   All variables names are laid out in the same manner as Forestry Canada 
  #   Fire Danger Group (FCFDG) (1992). Development and Structure of the 
  #   Canadian Forest Fire Behavior Prediction System." Technical Report 
  #   ST-X-3, Forestry Canada, Ottawa, Ontario.
  #
  # Args:
  #   FUELTYPE: The Fire Behaviour Prediction FuelType
  #   ROSeq:    The predicted equilibrium rate of spread (m/min)
  #   HR (t):   The elapsed time (min)
  #   CFB:      Crown Fraction Burned
  #   
  # Returns:
  #   DISTt:    Head fire spread distance at time t
  #
  #############################################################################
  #Eq. 72 (FCFDG 1992)
  #Calculate the alpha constant for the DISTt calculation
  alpha <- ifelse(FUELTYPE %in% c("C1", "O1A", "O1B", "S1", "S2", "S3", "D1"),
                  0.115,
                  0.115 - 18.8 * (CFB**2.5) * exp(-8* CFB))
  #Eq. 71 (FCFDG 1992) Calculate Head fire spread distance
  DISTt  <- ROSeq * (HR + exp(-alpha * HR) / alpha - 1 / alpha)
  
  return(DISTt)
}
