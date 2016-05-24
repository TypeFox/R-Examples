FitEigenValues <- function(rcov, phiGrid, phi, noEig) {
  buff <- .Machine$double.eps * max(abs(phiGrid)) * 3  
  
  if (is.null(noEig)) 
    noEig <- ncol(phi)
    
  # Get design matrix X:
  X <- apply(phi[, 1:noEig], 2, function(y) 
    approx(phiGrid, y, rcov$tPairs[, 1])$y * approx(phiGrid, y, rcov$tPairs[, 2])$y
  )
  
  if (class(rcov) == 'RawCov') {
    dat <- cbind(y=rcov[['cxxn']], X)
    dat <- dat[rcov$tPairs[, 1] > min(phiGrid) - buff & 
      rcov$tPairs[, 1] < max(phiGrid) + buff & 
      rcov$tPairs[, 2] > min(phiGrid) - buff & 
      rcov$tPairs[, 2] < max(phiGrid) + buff, ]
    mod <- lm(y ~ . - 1, data.frame(dat))
  } else if (class(rcov) == 'BinnedRawCov') {
    dat <- cbind(y=rcov[['meanVals']], X)
    dat <- dat[rcov$tPairs[, 1] > min(phiGrid) - buff & 
      rcov$tPairs[, 1] < max(phiGrid) + buff & 
      rcov$tPairs[, 2] > min(phiGrid) - buff & 
      rcov$tPairs[, 2] < max(phiGrid) + buff, ]    
    mod <- lm(y ~ . - 1, data.frame(dat), weights=rcov[['count']])
  }
  
  lam <- unname(mod[['coefficients']])
  
  if (any(lam <= 0))
    warning('Fit method produces negative estimates of eigenvalues')
    
  return(lam)
}
