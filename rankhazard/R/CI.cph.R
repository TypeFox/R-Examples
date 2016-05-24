rankhazard_CI.cph <- function(cphobj,  x, refpoints, CI_level){
  
  ### Calculates the values for the confidence intervals ###
  
  factorlevels = cphobj$Design$parms
  covariatelabs = cphobj$Design$name
  
  coefs <- cphobj$coef
  coefs_low <- confint(cphobj, level = CI_level)[, 1]
  coefs_upp <- confint(cphobj, level = CI_level)[, 2]
  
  Values <- rankhazard_CI(modelobj = cphobj, x = x, coefs = coefs, refpoints = refpoints, 
                                  factorlevels = factorlevels, covariatelabs = covariatelabs)
  CIlow <- rankhazard_CI(modelobj = cphobj, x = x, coefs = coefs_low, refpoints = refpoints, 
                                 factorlevels = factorlevels, covariatelabs = covariatelabs)
  CIupp <- rankhazard_CI(modelobj = cphobj, x = x, coefs = coefs_upp, refpoints = refpoints, 
                                 factorlevels = factorlevels, covariatelabs = covariatelabs)
  
  return(list(x = Values$x, xp = Values$xp, refvalues = Values$refvalues, low = CIlow$xp, lowrefvalues = CIlow$refvalues, upp = CIupp$xp, upprefvalues = CIupp$refvalues, select_CI = Values$select_CI))
  
}