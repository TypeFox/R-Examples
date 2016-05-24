#roxygen2 documentation in AceUnivariate.R
DeFriesFulkerMethod1 <- function( dataSet, oName_S1, oName_S2, rName="R" ) { 
  lmDetails <- stats::lm(
    dataSet[, oName_S1] ~ 
      1 + 
      dataSet[, oName_S2] +
      dataSet[, rName] + 
      dataSet[, oName_S2]*dataSet[, rName]
  )
  
  brief <- base::summary(lmDetails)
  coeficients <- stats::coef(brief)
  nDouble <- base::length(brief$residuals) 
  #b0 <- coeficients["(Intercept)", "Estimate"]
  b1 <- coeficients["dataSet[, oName_S2]", "Estimate"]  
  #b2 <- coeficients["R", "Estimate"]
  b3 <- coeficients["dataSet[, oName_S2]:dataSet[, rName]", "Estimate"]
  eSquared <- 1 - (b1+b3)
  
  details <- base::list(lm=lmDetails)
  aceEstimate <- NlsyLinks::CreateAceEstimate(aSquared=b3, cSquared=b1, eSquared=eSquared, caseCount=nDouble, details=details)
  return( aceEstimate )
}

#roxygen2 documentation in AceUnivariate.R
DeFriesFulkerMethod3 <- function( dataSet, oName_S1, oName_S2, rName="R" ) { 
  dv_S1Centered <- base::scale(dataSet[, oName_S1], center=TRUE, scale=FALSE)
  dv_S2Centered <- base::scale(dataSet[, oName_S2], center=TRUE, scale=FALSE)
  interaction <- dv_S2Centered*dataSet[, rName]
    
  lmDetails <- stats::lm(
    dv_S1Centered ~ 
      0 + 
      dv_S2Centered + 
      interaction
  ) #The '0' specifies and intercept-free model.
  brief <- base::summary(lmDetails)
    
  coeficients <- stats::coef(brief)
  nDouble <- base::length(brief$residuals) 
  b1 <- coeficients["dv_S2Centered", "Estimate"]  
  b2 <- coeficients["interaction", "Estimate"]
  eSquared <- 1 - (b1+b2)
    
  details <- base::list(lm=lmDetails)
  aceEstimate <- NlsyLinks::CreateAceEstimate(aSquared=b2, cSquared=b1, eSquared=eSquared, caseCount=nDouble, details=details)
  return( aceEstimate )
}
