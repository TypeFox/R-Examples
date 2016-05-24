.LBtcalc <- function(FUELTYPE, LB, HR, CFB) {
  #############################################################################
  # Description:
  #   Computes the Length to Breadth ratio of an elliptically shaped fire at
  #   elapsed time since ignition. Equations are from listed FCFDG (1992) and
  #   Wotton et. al. (2009), and are marked as such.
  #
  #   All variables names are laid out in the same manner as Forestry Canada 
  #   Fire Danger Group (FCFDG) (1992). Development and Structure of the 
  #   Canadian Forest Fire Behavior Prediction System." Technical Report 
  #   ST-X-3, Forestry Canada, Ottawa, Ontario.
  #
  #   Wotton, B.M., Alexander, M.E., Taylor, S.W. 2009. Updates and revisions to
  #   the 1992 Canadian forest fire behavior prediction system. Nat. Resour. 
  #   Can., Can. For. Serv., Great Lakes For. Cent., Sault Ste. Marie, Ontario, 
  #   Canada. Information Report GLC-X-10, 45p.
  #
  # Args:
  #   FUELTYPE: The Fire Behaviour Prediction FuelType
  #         LB: Length to Breadth ratio
  #         HR: Time since ignition (hours)
  #        CFB: Crown Fraction Burned
  # Returns:
  #   LBt: Length to Breadth ratio at time since ignition
  #
  #############################################################################
  #
  #Eq. 72 (FCFDG 1992) - alpha constant value, dependent on fuel type
  alpha <- ifelse(FUELTYPE %in% c("C1", "O1A", "O1B", "S1", "S2", "S3", "D1"),
              0.115,
              0.115 - 18.8 * (CFB**2.5) * exp(-8* CFB))  
  #Eq. 81 (Wotton et.al. 2009) - LB at time since ignition
	LBt <- (LB - 1) * (1 - exp(-alpha * HR)) + 1
  return(LBt)
}
