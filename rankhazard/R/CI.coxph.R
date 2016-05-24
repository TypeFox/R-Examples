rankhazard_CI.coxph <- function(coxphobj, x, refpoints, CI_level){
  
  ### Calculates the values for the confidence intervals ###
  
  factorlevels = coxphobj$xlevels
  covariatelabs = attr(coxphobj$terms, "term.labels")
  
  coefs <- coxphobj$coef
  coefs_low <- confint(coxphobj, level = CI_level)[, 1]
  coefs_upp <- confint(coxphobj, level = CI_level)[, 2]
  
  Values <- rankhazard_CI(modelobj = coxphobj, x = x, coefs = coefs, refpoints = refpoints, 
                                  factorlevels = factorlevels, covariatelabs = covariatelabs)
  CIlow <- rankhazard_CI(modelobj = coxphobj, x = x, coefs = coefs_low, refpoints = refpoints, 
                                 factorlevels = factorlevels, covariatelabs = covariatelabs)
  CIupp <- rankhazard_CI(modelobj = coxphobj, x = x, coefs = coefs_upp, refpoints = refpoints, 
                                 factorlevels = factorlevels, covariatelabs = covariatelabs)
  
  return(list(x = Values$x, xp = Values$xp, refvalues = Values$refvalues, low = CIlow$xp, lowrefvalues = CIlow$refvalues, upp = CIupp$xp, upprefvalues = CIupp$refvalues, select_CI = Values$select_CI))
  
}

