estimatePsi <- function(theta, theta0, fX, obs, sigma2, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ## Argument 'theta':
  mode <- mode(theta);
  if (mode != "function") {
    throw("Argument 'theta' should be of mode 'function', not '", mode);
  }

  ## Argument 'theta0':
  mode <- mode(theta0);
  if (mode != "function") {
    throw("Argument 'theta0' should be of mode 'function', not '", mode);
  }

  ## Argument 'fX':
  mode <- mode(fX);
  if (mode != "function") {
    throw("Argument 'fX' should be of mode 'function', not '", mode);
  }
  
  ## Argument 'obs':
  obs <- validateArgumentObs(obs);
  
  ## Argument 'sigma2':
  sigma2 <- Arguments$getNumeric(sigma2);


  T <- theta(obs[, c("X", "W")]);
  verbose && cat(verbose, "theta(X, W):");
  verbose && str(verbose, T);
  T0 <- theta0(obs[, "W", drop=FALSE]);
  verbose && cat(verbose, "theta0(W):");
  verbose && str(verbose, T0);

  mean.psi1 <- mean(fX(obs) * (T - T0))/sigma2;
  sd.psi1 <- sd(fX(obs) * (T - T0))/sigma2;
    
  list(mean=mean.psi1, sd=sd.psi1);
}


############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

