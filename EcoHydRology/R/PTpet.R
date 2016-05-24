PTpet <-
function(Rn, T_C, PTconstant = 1.26){
    LatentHtEvap <- 2500	# kJ/kg
    DensityWater <- 1000	# kg/m3
    PsychConstant <- 0.066	# kPa/K
    potentialET <- PTconstant * SatVapPresSlope(T_C) * Rn/((SatVapPresSlope(T_C) + PsychConstant) * (LatentHtEvap * DensityWater))
    potentialET[which(potentialET < 0)] <- 0
    return(signif(potentialET,2))

}
