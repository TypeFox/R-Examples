.TFCcalc <- function(FUELTYPE, CFL, CFB, SFC, PC, PDF, option = "TFC") {
  #############################################################################
  # Description:
  #   Computes the Total (Surface + Crown) Fuel Consumption by Fuel Type.
  #   All variables names are laid out in the same manner as FCFDG (1992) or
  #   Wotton et. al (2009) 
  
  #   Forestry Canada Fire Danger Group (FCFDG) (1992). "Development and 
  #   Structure of the Canadian Forest Fire Behavior Prediction System." 
  #   Technical Report ST-X-3, Forestry Canada, Ottawa, Ontario.
  #
  #   Wotton, B.M., Alexander, M.E., Taylor, S.W. 2009. Updates and revisions to
  #   the 1992 Canadian forest fire behavior prediction system. Nat. Resour. 
  #   Can., Can. For. Serv., Great Lakes For. Cent., Sault Ste. Marie, Ontario, 
  #   Canada. Information Report GLC-X-10, 45p.
  #
  # Args:
  #   FUELTYPE: The Fire Behaviour Prediction FuelType
  #        CFL: Crown Fuel Load (kg/m^2)
  #        CFB: Crown Fraction Burned (0-1)
  #        SFC: Surface Fuel Consumption (kg/m^2)
  #         PC: Percent Conifer (%)
  #        PDF: Percent Dead Balsam Fir (%)
  #     option: Type of output (TFC, CFC, default=TFC)
  # Returns:
  #        TFC: Total (Surface + Crown) Fuel Consumption (kg/m^2)
  #       OR
  #        CFC: Crown Fuel Consumption (kg/m^2)
  # #############################################################################
  #Eq. 66a (Wotton 2009) - Crown Fuel Consumption (CFC)
  CFC <- CFL * CFB
  #Eq. 66b (Wotton 2009) - CFC for M1/M2 types
  CFC <- ifelse(FUELTYPE %in% c("M1", "M2"), PC / 100 * CFC, CFC)
  #Eq. 66c (Wotton 2009) - CFC for M3/M4 types
  CFC <- ifelse(FUELTYPE %in% c("M3","M4"), PDF / 100 * CFC, CFC)
  #Return CFC if requested
  if (option=="CFC")
    return(CFC)
  #Eq. 67 (FCFDG 1992) - Total Fuel Consumption
  TFC <- SFC + CFC
  return(TFC)
}
