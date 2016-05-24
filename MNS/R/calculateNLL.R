calculateNLL <-
function(resp, des, fixEf, reVar, resVar){
  # Calculate the negative log-likelihood (NLL) of the current LME fit
  # Due to properties of EM algorithm this should always be DECREASING (note that M in our case in MINIMIZATION)
  #
  # INPUT:
  #      - resp: vector of responses (assumed for one subject i.e., X^{(i)}_{v})
  #      - des: fixed design for one subject (i.e., X^{(i)}_{\backslash v})
  #      - fixEf: regression coefficient for fixed effects
  #      - reVar: vector of random effect variances
  #      - resVar: scalar for residual standard deviation
  #
  #
  
  n = nrow(des)
  
  #V = (resVar**2) * (des %*% diag(reVar**2) %*% t(des) + diag(n))
  V =  (resVar**2) *(des %*% diag(reVar) %*% diag(reVar) %*% t(des) + diag(n))
  
  #logDet = 0.5*log(det(V))
  #Res = 0.5 * t(resp - des %*% fixEf) %*% solve(V) %*% (resp-des %*% fixEf)
  
  return(dmvnorm(as.vector(resp), mean = des %*% fixEf, sigma = V, log = TRUE))
  #return(Res)
  #return(as.numeric(logDet+Res))
}
