# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


PSTM = function (ffs, f0, template, winner = TRUE){
  ffs = as.matrix(ffs)
  if (length(ffs) == 3 & ncol(ffs) == 1) ffs = t(ffs)
  if (ncol(ffs) < 3) stop ('At least three formants must be provided.')
  n = nrow (ffs)
  f0 = as.matrix(f0)
  if (nrow(ffs) != length(f0)) stop ('Rows in ffs must equal length of f0.')
    
  if (missing(template)){
    covariance = matrix (c(1.471405e-02,1.042306e-03, -1.011399e-03,  1.042306e-03,
    1.005482e-02,  5.882612e-06, -1.011399e-03,  5.882612e-06,  5.674205e-03), 3, 3)
    
    means =  matrix (c(-1.4191406,0.62950438,0.8487829,-0.9816006,0.41787938,
    0.7149399,-0.9211756,0.49334438,0.7382749,-0.7035106,0.32145438,0.7001354,-0.4457406,0.22904938,
    0.6735179,-0.5450806,-0.09682062,0.6557979,-0.6572806,0.01781438,0.6558074,-0.8547056,-0.17574562,
    0.6554719,-0.9187456,-0.09007562,0.6727639,-1.2477756,0.01196438,0.6208954), 10, 3, byrow = TRUE)
    
	rownames (means) = c('i','I','e','E','{','A','V','o','U','u')
	
    template = list (means = means, covariance = covariance)
    class (template) = 'template'
  }
  if (class(template) != 'template') stop ('Innapropriate template provided.')

  f0weight=2; covarweight=1; psiprioweight=1;  

  means = as.matrix(template$means)
  if (ncol(ffs) != ncol(means)) stop ('Same number of formants must be provided in ffs and means.')
  covariance = as.matrix(template$covariance)
  if (ncol(covariance) != ncol(means)) stop ('Same number of formants must be provided in covariance and means.')
   
  if (mean(ffs) > 20) ffs = log(ffs)
  if (mean(f0) > 20) f0 = log(f0)
  
  psiPriorMean = 7.2333
  psiPriorVar = (.1284)^2 * 1
  
  f0PsiSlope = 2.14452
  f0PsiIntercept = - 10.3233
  f0VarGivenPsi = (.132703)^2
  
  covariance = solve (covariance * covarweight)
  if (winner == TRUE) output = data.frame()
  if (winner == FALSE) output = list()
  for (i in 1:n){
    outpsi = rep (0, nrow(means))
    dvowels = rep (0, nrow(means))
    ## merged derivative of log-likelihood lines
    ## slope (sd) is constant, only intercept changes
    slope = -sum(covariance) + (-f0PsiSlope/f0VarGivenPsi) + (-1/psiPriorVar)
    baseintercept = f0[i]/f0VarGivenPsi - (f0PsiIntercept/f0VarGivenPsi) + (psiPriorMean/psiPriorVar)
    k = ncol (ffs)  ## this is the multivariate normal density given means and an observation
  
    for (j in 1:nrow(means)){
      intercept = baseintercept + sum ((ffs[i,]-means[j,])%*%covariance) 
      outpsi[j] = -intercept/slope
      normd = ffs[i,] - outpsi[j]
      dvowels[j] = exp( -.5 * ( (normd - means[j,]) %*% covariance %*% (normd - means[j,]) ) )
    }
    dvowels = dvowels / sum (dvowels)
    tmp = data.frame (vowel = rownames(means), psi = outpsi, postprob = dvowels)
    
    if (winner == FALSE) output[[i]] = tmp
    if (winner == TRUE) output = rbind (output, tmp[order(tmp[,3], decreasing = TRUE)[1],])
  }
  if (winner == TRUE) rownames (output) = 1:nrow(output)
  return (output)
}
