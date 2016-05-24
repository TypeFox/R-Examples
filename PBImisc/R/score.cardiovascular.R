#age - age
#cholesterol - in mmol/L
#SBP - in mmHg
#currentSmoker - 1 for current smokers, 0 fof non smokers
  

calculateRisk <- function(age, cholesterol, SBP, currentSmoker, betaSmoker, betaSBP, betaChol, coefs) {
  # step 1 risks
  Sage0 = exp(-exp(coefs["alpha"])*(age - 20)^coefs["p"])
  Sage10 = exp(-exp(coefs["alpha"])*(age - 10)^coefs["p"])
  # step 2 weights
  w = betaChol*(cholesterol - 6) + betaSBP*(SBP - 120) + betaSmoker
  # step 3 weighted risks
  Sage   = (Sage0)^exp(coefs["w"]) 
  Sage1 = (Sage10)^exp(coefs["w"]) 
  # step 4 - 10 years survival
  S10 = Sage1/Sage
  # step 5 - endpoint
  Risk10 = 1 - S10
  Risk10
}

calculateScoreEur <- function(age, cholesterol, SBP, currentSmoker, gender = "Men", risk = "Low risk") {
  betaSmoker = c(0.71, 0.63)
  betaSBP    = c(0.018, 0.022)
  betaChol   = c(0.24, 0.02)
  
  coeffs <- array(c(-22.1, 4.71, -26.7, 5.64, -29.8, 6.36, -31.0, 6.62, -21.0, 4.62, -25.7, 5.47, -28.7, 6.23, -30.0, 6.42), 
                  c(2,2,2,2),
                  dimnames = list(c("alpha", "p"), c("CHD", "non CHD"), c("Men", "Women"), c("Low risk", "High risk")))
  
  # step 6 - score
  CVDrisk = calculateRisk(age, cholesterol, SBP, currentSmoker,
                          betaSmoker[1], betaSBP[1], betaChol[1], coeffs[,"CHD",gender,risk])
  NonCVDrisk = calculateRisk(age, cholesterol, SBP, currentSmoker,
                          betaSmoker[2], betaSBP[2], betaChol[2], coeffs[,"non CHD",gender,risk])
 
  CVDrisk + NonCVDrisk
}
