##=======================================================================
## auxiliary variables for a gumbel distribution
## see "Renext Computing Details".
##
##=======================================================================

gumbelAux <- function(r) {
  res1 <- rep(0, length(r))
  res2 <- rep(pi * pi /6, length(r))
  for (i in 1:length(r)) {
    if (r[i] > 1) {
      temp <- 1 / (1:(r[i]-1)) 
      res1[i] <- sum(temp)
      res2[i] <- res2[i] - sum(temp^2) 
    }
  }
  list(r = r, nu = res1, kappa2 = res2)
}

##=======================================================================
## compute initial (cheap estimate) for 'lambda' and 'sigma'
## for the aggregated POT model with exponential exceedances
## using MAX historical data only.
## 
## Here 'MAX' is a list as computed with the 'makeMAXdata' function
##
##============================================================1===========

parIni.MAX <- function(MAX, threshold, distname.y = "exp") {

  if (!(distname.y %in% c("exp", "exponential", "gpd", "GPD"))) {
    stop("unauthorised value for 'distname.y'")
  }
  r <- MAX$r
  nBloc <- length(r)
  if (nBloc < 2L) stop("The number of blocks must be >= 2") 
  
  ## minimal r-largest statistic
  y <- unlist(lapply(MAX$data, min))
  Spac <- unlist(lapply(MAX$data, spacings))
  nSpac <- length(Spac) 
  aux <- gumbelAux(r)
  
  ## compute sigmaHat
  sigmaSpac <- mean(Spac)
  
  ## regression estimates, see "Renext Computing Details"
  x <- log(MAX$effDuration) - aux$nu
  wt <-  1.0 / aux$kappa2
  covmat <- cov.wt(cbind(x, y), wt = wt, method = "ML")$cov
  
  a <- 1.0 / mean(wt)
  betaReg <- polyroot(c(-covmat[2L, 2L], covmat[1L, 2L], a))
  betaReg <- sort(Re(betaReg))[2L]

  if (nSpac) {
    den <- nSpac + nBloc / 2.0
    wSpac <- nSpac / den
    sigmaHat <- wSpac * sigmaSpac  + (1 - wSpac) * betaReg
  } else {
    sigmaHat <- betaReg
  }
    
  alphaReg <- weighted.mean(y, w = wt) - betaReg * weighted.mean(x, w = wt)
  lambdaHat <- exp( (alphaReg + digamma(1) * sigmaHat - threshold) / sigmaHat )

 if (distname.y %in% c("exp", "exponential")) {
    return(c("lambda" = lambdaHat, "rate" = 1.0 / sigmaHat))
    ## "alphaReg" = alphaReg, "betaReg" = betaReg, "sigmaSpac" = sigmaSpac)
  } else if (distname.y %in% c("gpd", "GPD")) {
    return(c("lambda" = lambdaHat, "scale" = sigmaHat, "shape" = 0.0))
  }
 
  
}

##=======================================================================
## compute initial (cheap estimate) for 'lambda' and 'sigma'
## for the aggregated POT model with exponential exceedances
## using OTS data only.
## 
## Here 'OTS' is a list as computed with the 'makeOTSdata' function
##
##=======================================================================

parIni.OTS <- function(OTS, threshold, distname.y = "exp") {

  if (!(distname.y %in% c("exp", "exponential", "gpd"))) {
    stop("unauthorised value for 'distname.y'")
  }

  dThreshold <- OTS$threshold - threshold
    
  if (any(dThreshold <= 0)) {
    stop("all values in 'OTSthreshold' must be > 'threshold'")
  }
  
  ## estimate the rate using a GLM, see "Renext Computing Details"
  R <-  OTS$r 
  lw <- log(OTS$effDuration)
  fit <- glm(R ~ dThreshold, family = poisson, offset = lw)
  betaReg <- as.numeric(coef(fit))
  
  lambdaReg <- exp(betaReg[1])
  ## sigmaReg <- - 1.0 / betaReg[2] ## not used 
  
  Ybar <- unlist(lapply(OTS$data, mean)) - unlist(OTS$threshold)
  sigmaHat <- weighted.mean(Ybar, w = OTS$r)

  if (distname.y %in% c("exp", "exponential")) {
    return(c("lambda" = lambdaReg, "rate" = 1.0 / sigmaHat)) 
  } else if (distname.y == "gpd") {
    return(c("lambda" = lambdaReg, "scale" = sigmaHat, "shape" = 0.0))
  }
  
}


