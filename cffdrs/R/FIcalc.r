.FIcalc <- function(FC, ROS){
  #############################################################################
  # Description:
  #   Calculate the Predicted Fire Intensity
  #
  #   All variables names are laid out in the same manner as Forestry Canada 
  #   Fire Danger Group (FCFDG) (1992). Development and Structure of the 
  #   Canadian Forest Fire Behavior Prediction System." Technical Report 
  #   ST-X-3, Forestry Canada, Ottawa, Ontario.
  #
  # Args:
  #   FC:   Fuel Consumption (kg/m^2)
  #   ROS:  Rate of Spread (m/min)
  #   
  # Returns:
  #   FI:   Fire Intensity (kW/m)
  #
  #############################################################################
  #Eq. 69 (FCFDG 1992) Fire Intensity (kW/m)
  FI   <- 300 * FC * ROS
  return(FI)
}
