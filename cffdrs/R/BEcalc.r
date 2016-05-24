.BEcalc<-function(FUELTYPE, BUI) {
  #############################################################################
  # Description:
  #   Computes the Buildup Effect on Fire Spread Rate. 
  #
  #   All variables names are laid out in the same manner as Forestry Canada 
  #   Fire Danger Group (FCFDG) (1992). Development and Structure of the 
  #   Canadian Forest Fire Behavior Prediction System." Technical Report 
  #   ST-X-3, Forestry Canada, Ottawa, Ontario.
  #
  # Args:
  #   FUELTYPE: The Fire Behaviour Prediction FuelType
  #   BUI:      The Buildup Index value
  # Returns:
  #   BE: The Buildup Effect
  #
  #############################################################################
  #Fuel Type String represenations
  d <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "D1", "M1", "M2", "M3",
         "M4","S1", "S2", "S3", "O1A", "O1B")
  #The average BUI for the fuel type - as referenced by the "d" list above
  BUIo <- c(72, 64, 62, 66, 56, 62, 106, 32, 50, 50, 50, 50, 38, 63, 31, 01, 
            01)
  #Proportion of maximum possible spread rate that is reached at a standard BUI
  Q <- c(0.9, 0.7, 0.75, 0.8, 0.8, 0.8, 0.85, 0.9, 0.8, 0.8, 0.8, 0.8, 0.75, 
         0.75, 0.75, 1.0, 1.0)
  names(BUIo) <- names(Q)<-d
  
  #Eq. 54 (FCFDG 1992) The Buildup Effect
  BE<- ifelse(BUI > 0 & BUIo[FUELTYPE] > 0,
    exp(50 * log(Q[FUELTYPE]) * (1 / BUI - 1 / BUIo[FUELTYPE])), 1)
  
  return(as.numeric(BE))
}