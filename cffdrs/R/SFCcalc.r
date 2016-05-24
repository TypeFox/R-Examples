.SFCcalc <- function(FUELTYPE, FFMC, BUI, PC, GFL) {
  #############################################################################
  # Description:
  #   Computes the Surface Fuel Consumption by Fuel Type.
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
  #        BUI: Buildup Index
  #       FFMC: Fine Fuel Moisture Code
  #         PC: Percent Conifer (%)
  #        GFL: Grass Fuel Load (kg/m^2)
  # Returns:
  #        SFC: Surface Fuel Consumption (kg/m^2)
  #
  #############################################################################
  SFC <- rep(-999,length(FFMC))
  #Eqs. 9a, 9b (Wotton et. al. 2009) - Solving the lower bound of FFMC value
  # for the C1 fuel type SFC calculation
  SFC <- ifelse(FUELTYPE=="C1", 
          ifelse(FFMC > 84, 
            0.75 + 0.75 * (1 - exp(-0.23 * (FFMC - 84)))**0.5,
            0.75 - 0.75 * (1 - exp(-0.23 * (84 - FFMC)))**0.5),
          SFC)
  #Eq. 10 (FCFDG 1992) - C2, M3, and M4 Fuel Types
  SFC <- ifelse(FUELTYPE == "C2" | FUELTYPE == "M3" | FUELTYPE == "M4", 
          5.0 * (1 - exp(-0.0115 * BUI)),
          SFC)
  #Eq. 11 (FCFDG 1992) - C3, C4 Fuel Types
  SFC <- ifelse(FUELTYPE == "C3" | FUELTYPE == "C4", 
          5.0 * (1 - exp(-0.0164 * BUI))**2.24,
          SFC)
  #Eq. 12 (FCFDG 1992) - C5, C6 Fuel Types
  SFC <- ifelse(FUELTYPE == "C5" | FUELTYPE == "C6",
          5.0 * (1 - exp(-0.0149 * BUI))**2.48,
          SFC)
  #Eqs. 13, 14, 15 (FCFDG 1992) - C7 Fuel Types
  SFC <- ifelse(FUELTYPE == "C7", 
          ifelse(FFMC > 70, 2 * (1 - exp(-0.104 * (FFMC - 70))), 
            0) + 1.5 * (1 - exp(-0.0201 * BUI)),
          SFC)
  #Eq. 16 (FCFDG 1992) - D1 Fuel Type
  SFC <- ifelse(FUELTYPE == "D1", 1.5 * (1 - exp(-0.0183 * BUI)), SFC)
  #Eq. 17 (FCFDG 1992) - M1 and M2 Fuel Types
  SFC <- ifelse(FUELTYPE == "M1" | FUELTYPE == "M2", 
          PC / 100 * (5.0 * (1 - exp(-0.0115 * BUI))) + 
           ((100 - PC) / 100 * (1.5 * (1 - exp(-0.0183 * BUI)))), 
          SFC)
  #Eq. 18 (FCFDG 1992) - Grass Fuel Types
  SFC <- ifelse(FUELTYPE == "O1A" | FUELTYPE == "O1B", GFL, SFC)
  #Eq. 19, 20, 25 (FCFDG 1992) - S1 Fuel Type
  SFC <- ifelse(FUELTYPE == "S1", 
          4.0 * (1 - exp(-0.025 * BUI)) + 4.0 * (1 - exp(-0.034 * BUI)), 
          SFC)
  #Eq. 21, 22, 25 (FCFDG 1992) - S2 Fuel Type
  SFC <- ifelse(FUELTYPE == "S2", 
          10.0 * (1 - exp(-0.013 * BUI)) + 6.0 * (1 - exp(-0.060 * BUI)), 
          SFC)
  #Eq. 23, 24, 25 (FCFDG 1992) - S3 Fuel Type
  SFC <- ifelse(FUELTYPE == "S3", 
          12.0 * (1 - exp(-0.0166 * BUI)) + 20.0 * (1-exp(-0.0210 * BUI)),
          SFC)
  #Constrain SFC value
  SFC <- ifelse(SFC <= 0, 0.000001, SFC)
  
  return(SFC)
}
