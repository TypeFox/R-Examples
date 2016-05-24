GetRho <- function(y, t, optns, mu, obsGrid, fittedCov, lambda, phi, sigma2) {
  
  optnsTmp <- optns
  optnsTmp$verbose <- FALSE 
  for (j in 1:2) {
    yhat <- GetCEScores(y, t, optnsTmp, mu, obsGrid, fittedCov, lambda, phi, sigma2)[3, ] 
    sigma2 <- mean(mapply(function(a, b) mean((a - b)^2, na.rm=TRUE), yhat, y), na.rm=TRUE)
  }
    
  R <- sqrt((trapzRcpp(obsGrid, mu ^ 2) + sum(lambda)) / diff(range(obsGrid)))
  a1 <- 0.01; a2 <- 0.22
  etaCand <- seq(a1, a2, length.out=50)
  rhoCand <- etaCand * R
  rhoCand <- rhoCand[rhoCand > sigma2]
  rhoCand <- c(sigma2, rhoCand)
  
  leaveOutInd <- RandTime(t, isRandom=FALSE)

  cvScores <- sapply(rhoCand, cvRho, leaveOutInd=leaveOutInd, y=y, t=t, optns=optns, mu=mu, obsGrid=obsGrid, fittedCov=fittedCov, lambda=lambda, phi=phi)

  return(rhoCand[which.min(cvScores)])
}


# sigma2* = max(rho, sigma2) = rho
# Get the CV score for a given rho. Correspond to getScores2.m
cvRho <- function(rho, leaveOutInd, y, t, optns, mu, obsGrid, fittedCov, lambda, phi) {

  Sigma_Y <- fittedCov + diag(rho, nrow(phi))

  MuPhiSig <- GetMuPhiSig(t, obsGrid, mu, phi, Sigma_Y)

  yhat <- mapply(function(yVec, muphisig, ind) 
         GetIndCEScores(yVec, muphisig$muVec, lambda, muphisig$phiMat,
                        muphisig$Sigma_Yi, newyInd=ind)$fittedY, 
                 y, MuPhiSig, leaveOutInd) 
  
  yobs <- mapply(`[`, y, leaveOutInd)

  return(sum((na.omit(unlist(yobs)) - unlist(yhat))^2, na.rm=TRUE))
}


# sample one observation from each tVec in t. The 'non-random' sampling is for testing against Matlab.
RandTime <- function(t, isRandom=TRUE) {
  ni <- sapply(t, length)
  if (all(ni <= 1)) 
    stop('None of the individuals have >= 2 observations. Cannot use rho')

  if (isRandom) {
    ind <- sapply(ni, sample, size=1)
  } else {
    ind <- ((1000 + 1:length(t)) %% ni) + 1
  }

  return(ind)
}
