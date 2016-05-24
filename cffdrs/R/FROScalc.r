.FROScalc <- function(ROS, BROS, LB){
  #############################################################################
  # Description:
  #   Calculate the Flank Fire Spread Rate. 
  #
  #   All variables names are laid out in the same manner as Forestry Canada 
  #   Fire Danger Group (FCFDG) (1992). Development and Structure of the 
  #   Canadian Forest Fire Behavior Prediction System." Technical Report 
  #   ST-X-3, Forestry Canada, Ottawa, Ontario.
  #
  # Args:
  #   ROS:    Fire Rate of Spread (m/min)
  #   BROS:   Back Fire Rate of Spread (m/min)
  #   LB:     Length to breadth ratio
  #   
  
  # Returns:
  #   FROS:   Flank Fire Spread Rate (m/min)
  #
  #############################################################################
  #Eq. 89 (FCFDG 1992)
  FROS <- (ROS + BROS) / LB / 2
  return(FROS)
}
