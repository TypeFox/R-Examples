setMethodS3("updatePsi", "NPVI", function(this, B, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ## Argument 'B':
  B <- Arguments$getInteger(B);
  
  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose);

  ## Retrieve parameters
  family <- getFamily(this)
  fX <- getFX(this)
  obs <- getObs(this)
  W <- obs[, "W"]
  X <- fX(obs)
  Xq <- getXq(this)
  
  g <- getG(this);
  mu <- getMu(this);
  muAux <- getMuAux(this);
  sigma2 <- getSigma2(this);
  weightsW <- getWeightsW(this)

  ## Perform 'B' simulations according to the estimated parameters
  verbose && enter(verbose, "Simulating ", B, " observations");
  obsB <- simulateData(B, W, X, Xq, g, mu, muAux, sigma2, weightsW=weightsW, family=family, verbose=verbose)
  verbose && str(verbose, obsB);
  verbose && exit(verbose);

  ## Calculate 'theta' and 'theta0' on these B samples
  theta <- getTheta(this)
  theta0 <- getTheta0(this)

  ## Estimate psiPn:
  psi1 <- estimatePsi(theta=theta, theta0=theta0, fX=fX, obs=obs, sigma2=sigma2, verbose=verbose) 
  this$.psiPn <- psi1$mean;
  this$.psiPn.sd <- psi1$sd;
  ## Estimate psi:
  psi0 <- estimatePsi(theta=theta, theta0=theta0, fX=fX, obs=obsB, sigma2=sigma2, verbose=verbose) 
  this$.psi <- psi0$mean;
  this$.psi.sd <- psi0$sd;
  rm(obsB)
})

############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

