#==============================================================================
# File: fitModel.R
#
# Author: Nathan Morris
#         Yeunjoo Song
#
# Notes: Main driver of fitting the model.
#
# History: Initial implementation
#          Revision - yes Jan 2013
#==============================================================================

#------------------------------------------------------------------------------
# Fitting the Model.
#------------------------------------------------------------------------------
.fitModel = function(model, y, x, vc, step1OptimControl, startValueControl, step2OptimControl)
{  
  numTrait = ncol(y$yAll[[1]])
  numCov   = ncol(x$xAll[[1]])
  numVC    = length(vc$vcAll[[1]])
  
  numC = numCov * numTrait
  numV = numVC * (numTrait*(numTrait+1)/2)

  df0 = (numC+numV) - length(paramNames(model))

  # error - Not identifiable, more parameters in model than in saturated model!
  if( df0 < 0 )
  {
    myFittedModel = new("strumFittedModel",
                         myStrumModel  = model,
                         modelValidity = 1)

    return(myFittedModel)
  }

  # 1. Fit stage1 
  #---------------
  s1 = .fitModelStep1(y, x, vc, step1OptimControl)

  # 2. Fit stage2 
  #---------------
  s2 = .fitModelStep2(s1, model, startValueControl, step2OptimControl)

  # error - Parameters are not locally identifiable!
  if( df0 != ncol(s2$dfOrthog) )
  {
    myFittedModel = new("strumFittedModel",
                         myStrumModel  = model,
                         modelValidity = 2)

    return(myFittedModel)
  }

  # 3. Test the model
  #-------------------
  tOut = .testModel(model, s1, s2)

  fitVal = 0

  # Same number of parameters in model and saturated model.  No fit test!
  if( !is.data.frame(tOut$parEst) )
  {
    fitVal = 4
  } else if( !is.data.frame(tOut$fitChi) )
  {
    fitVal = 3
  }

  myFittedModel = new("strumFittedModel",
                       myStrumModel              = model,
                       modelValidity             = fitVal,
                       fittedParameters          = tOut$parEst,
                       fittedParametersCovMatrix = s2$thetaCov,
                       deltaParameters           = s1$delta,
                       deltaParametersCovMatrix  = s1$deltaCov,
                       parDiff                   = s2$diff,
                       parDiffCovMatrix          = s2$diffCov,
                       chiTestOut                = tOut$fitChi,
                       fitIndices                = tOut$fitIndx)

  return(myFittedModel)
}

#------------------------------------------------------------------------------
# Fitting the Model, stage 1
# - Form a limited information estimate for the saturated model.
#------------------------------------------------------------------------------
.fitModelStep1 = function(y, x, vc, step1OptimControl)
{
  CV = .estimateDeltaParameter(y, x, vc, step1OptimControl)
  De = .estimateDeltaDerivative(y, x, vc, CV$C, CV$V)
  JF = .estimateDeltaCovariance(De)

  delta = .Call("getDelta", t(CV$C), CV$V)

  .printInfoLine("Fitting model step 1", "Done", 50, 2)

  return(list(CV=CV, K=JF$K, delta=delta, deltaCov=JF$deltaCov, w=JF$w, Wi=De$W))
}

#------------------------------------------------------------------------------
# Fitting the Model, stage 2
# - Estimate parameters and asymptotic covariance matrix of the model.
#------------------------------------------------------------------------------
.fitModelStep2 = function(s1, model, startValueControl, step2OptimControl)
{
  startTheta = .generateStartValue(model, s1$delta, s1$w, startValueControl)
  optimTheta = .estimateThetaParameter(model, s1$delta, s1$w, startTheta, step2OptimControl)
  parCov     = .estimateThetaCovariance(model, s1, optimTheta)

  .printInfoLine("Fitting model step 2", "Done", 50, 2)

  return(list(theta    = optimTheta,
              thetaCov = parCov$thetaCov,
              diff     = parCov$diff,
              diffCov  = parCov$diffCov,
              dfOrthog = parCov$dfOrthog,
              U        = parCov$U))
}

#------------------------------------------------------------------------------
# Test the model fit.
#------------------------------------------------------------------------------
.testModel = function(model, s1, s2)
{
  # Calculate SE, CI, and p values of theta.
  #------------------------------------------------------------------------------
  theta    = s2$theta
  thetaCov = s2$thetaCov

  myeps = 1.0e-11 #.Machine$double.eps (2.220446e-16) or 1.0e-11

  thetaVar = diag(thetaCov)
  thetaVar[which(abs(thetaVar) < myeps)] = 0

  if( length(which(thetaVar < 0)) > 0 )
  {
    return(list(parEst  = NA,
                fitChi  = NA,
                fitIndx = NA))
  }
    
  standardError = sqrt(thetaVar)

  lowerCI = theta + qnorm(.025) * standardError
  upperCI = theta - qnorm(.025) * standardError

  z      = theta/standardError
  pValue = pnorm(-abs(z))

  positive = grep("<", paramNames(model))
  lengthNV = length(paramNames(model)) - length(positive)
  pValue[1:lengthNV] = pValue[1:lengthNV] * 2.0

  dfPar = data.frame(estimate = theta,
                     stdError = standardError,
                     lowerCI  = lowerCI,
                     upperCI  = upperCI,
                     pValue   = pValue)

  if( length(s1$delta) == length(theta) )
  {
    return(list(parEst  = dfPar,
                fitChi  = NA,
                fitIndx = NA))
  }

  # Fit check
  #-----------
  N       = s1$K # or max(s1$usedN)
  dfDelta = length(s1$delta)
  dfTheta = length(theta)
  
  # 1. Naive chi-square statistics of model
  #-----------------------------------------
  T1naive  = sum(s2$diff*s2$diff*s1$w) * (N-1)
  df1      = dfDelta - dfTheta

  # 2. Mean scaled chi-square statistics
  #--------------------------------------
  dfOrthog = s2$dfOrthog
  gamma    = s1$deltaCov
  Ugamma   = s2$U %*% gamma
  trUgamma = sum(diag(Ugamma))

  temp     = (t(dfOrthog) %*% gamma %*% dfOrthog)
  out      = qr(temp)
  rank     = out$rank

  kappa1 = trUgamma / rank
  T1mean = T1naive / kappa1

  # 3. Mean/variance scaled chi-square statistics
  #-----------------------------------------------
  dstar      = trUgamma^2 / sum(diag(Ugamma^2))
  df1meanvar = round(dstar)

  kappa2    = trUgamma / df1meanvar
  T1meanvar = T1naive / kappa2

  pVal1naive   = pchisq(T1naive,   df1,        lower.tail = FALSE)
  pVal1mean    = pchisq(T1mean,    rank,       lower.tail = FALSE)
  pVal1meanvar = pchisq(T1meanvar, df1meanvar, lower.tail = FALSE)

  # 3. Corrected chi-square statistics
  #------------------------------------
  chiCorrect   = .getTheoCorrectT(s2$diffCov, T1naive, df1, s1$w, myeps)
  pVal1correct = chiCorrect$pVal
  T1correct    = chiCorrect$Tcorrect
  
  # Fit Test Results
  #------------------
  dfTest = data.frame(kappa   = c(1,          kappa1,    kappa2,       NA),
                      chiStat = c(T1naive,    T1mean,    T1meanvar,    T1correct),
                      df      = c(df1,        rank,      df1meanvar,   df1),
                      pValue  = c(pVal1naive, pVal1mean, pVal1meanvar, pVal1correct))
  
  row.names(dfTest) = c("Un-adjusted", "Mean adjusted", "Mean-Variance adjusted", "Theoretically corrected")

  # Additional Fit Indices
  #------------------------
  
  # 4. Comparative Fit Index (CFI)
  #--------------------------------
  chi0  = .getChiNullCorrect(s1, myeps)
  
  dNull = chi0$T0correct - chi0$df0
  dAlt  = T1correct - df1
  
  cfi = (dNull-dAlt) / dNull

  if( is.na(cfi) )     cfi = 1.0
  else if( cfi > 1.0 ) cfi = 1.0
  else if( cfi < 0.0 ) cfi = 0.0
 
  # 5. Root Mean Square Error Approximation (RMSEA)
  # - more theoretic work needed since we need to find
  #   a way to calculate the effective sample size N
  #-------------------------------------------------
  #rmsea = sqrt((T1correct-df1)/(df1*effectiveN))
  
  fitIndx = data.frame(value = c(cfi))
  row.names(fitIndx) = c("Comparative Fit Index (CFI)")

  .printInfoLine("Testing model fit", "Done", 50, 2)

  return(list(parEst  = dfPar,
              fitChi  = dfTest,
              fitIndx = fitIndx))
}

#------------------------------------------------------------------------------
# A function to calculate the theoretically corrected chi-square T
#  for fit test 
#------------------------------------------------------------------------------
.getTheoCorrectT = function(M, Tnaive, dfnaive, w, myeps)
{
  eM = eigen(M, symmetric=TRUE)
  eM$values[which(abs(eM$values) < myeps)] = 0
  eM$values[which(eM$values < 0)] = 0 
  eM$vectors[which(abs(eM$vectors) < myeps)] = 0
  
  eMsqrt = eM$vectors %*% diag(sqrt(eM$values)) %*% t(eM$vectors)
  
  B = t(eMsqrt) %*% diag(w) %*% eMsqrt
  
  eB = eigen(B, symmetric=TRUE)
  
  nEival = length(eB$values)
  nsims  = 10000
  
  simVals = matrix(rchisq(nsims * nEival, df=1), nrow=nsims, ncol=nEival)
  simVals = simVals %*% eB$values
  
  #pVal = mean(simVals < T0)
  pVal = mean(simVals > Tnaive)
  #if( pVal < (50/nsims) )
  if( pVal == 0 )
    Tcorrect = NA
  else
    Tcorrect = qchisq(pVal, dfnaive, lower.tail = FALSE)
  
  return(list(Tcorrect=Tcorrect, pVal=pVal))
}

#------------------------------------------------------------------------------
# A function which calculates the vector of partial derivatives of a function
#  numerically 
#------------------------------------------------------------------------------
.getChiNullCorrect = function(s1, myeps)
{
  # Find the positon of the parameters in null model
  # - null model only include intercept and variance
  C0  = matrix(0, nrow = nrow(s1$CV$C), ncol = ncol(s1$CV$C), byrow = FALSE)
  C0[1,] = seq(from=1, to=ncol(C0))
  V0 = list()  	
  for( j in 1:length(s1$CV$V) )
  {
    V0[[j]] = matrix(0, nrow(s1$CV$V[[1]]), ncol(s1$CV$V[[1]]))
    diag(V0[[j]]) = seq(from=(ncol(C0)*j+1), to=(ncol(C0)*j+ncol(C0)))
  }
  
  delta0pos = .Call("getDelta", t(C0), V0)
  
  dFdTheta0 = sapply(delta0pos, function(i) {as.numeric(1:max(delta0pos)==i)})
  
  # Calculate covariance matrix of null theta.
  #--------------------------------------------
  W = diag(s1$w)
  WJFJW = W %*% s1$deltaCov %*% W 
  
  A = dFdTheta0 %*% W     %*% t(dFdTheta0)
  B = dFdTheta0 %*% WJFJW %*% t(dFdTheta0)
  
  Ai = tryCatch(solve(A), error=function(e) return(ginv(A)))
  
  theta0Cov = (Ai %*% B %*% Ai) / s1$K
  
  # Calculate covariance matrix of diff=(delta-theta0).
  #--------------------------------------------------------
  diff0    = s1$delta - s1$delta*(delta0pos>0)
  dFAidFW  = t(dFdTheta0) %*% Ai %*% dFdTheta0 %*% W
  IdFAidFW = diag(nrow(W)) - dFAidFW
  
  diff0Cov = (IdFAidFW %*% s1$deltaCov %*% t(IdFAidFW))

  # Naive chi-square statistics of null model
  #-----------------------------------------
  T0naive = sum(diff0*diff0*s1$w) * (s1$K-1)
  df0     = length(s1$delta) - max(delta0pos)
  
  # Corrected chi-square statistics of null model
  #-----------------------------------------
  chiCorrect = .getTheoCorrectT(diff0Cov, T0naive, df0, s1$w, myeps)
  
  return(list(T0correct=chiCorrect$Tcorrect, df0=df0))
}
