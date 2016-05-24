stddev.fd <- function(fdobj)
{
  #  Compute the standard deviation functions for functional observations
  #  Argument:
  #  fdobj    ... a functional data object
  #  Return:
  #  STDFD ... a functional data object for the standard deviation functions

  #  Last modified 20081004;  previously modified 26 February 2007

  if (!(inherits(fdobj, "fd"))) stop(
		"Argument  fdobj not a functional data object.")

  coef     <- fdobj$coefs
  coefd    <- dim(coef)
  ndim     <- length(coefd)
  if (coefd[1] == 1) stop("Only one replication found.")
  nrep     <- coefd[2]
  ones     <- rep(1,nrep)

  basisobj <- fdobj$basis
  fdnames  <- fdobj$fdnames
  nbasis   <- basisobj$nbasis
  rangeval <- basisobj$rangeval

  neval   <- max(c(201,10*nbasis + 1))
  evalarg <- seq(rangeval[1], rangeval[2], length=neval)
  fdarray <- eval.fd(evalarg, fdobj)

  if (ndim == 2) {
    mnvec  <- (fdarray %*% ones)/nrep
    resmat <- fdarray - (c(mnvec) %o% ones)
    varvec <- (resmat^2 %*% ones)/(nrep-1)
    stdmat <- sqrt(varvec)
  } else {
    nvar <- coefd[3]
    stdmat <- matrix(0, neval, nvar)
    for (j in 1:nvar) {
      mnvecj  <- (fdarray[,,j] %*% ones)/nrep
      resmatj <- fdarray[,,j] - outer(c(mnvecj), ones)
      varvecj <- (resmatj^2 %*% ones)/(nrep-1)
      stdmatj <- sqrt(varvecj)
      stdmat[,j] <- stdmatj
    }
  }
  stdcoef <- project.basis(stdmat, evalarg, basisobj)
  names(fdnames)[2] <- "Std. Dev."
  names(fdnames)[3] <- paste("Std. Dev.",names(fdnames)[3])
  stdfd <- fd(stdcoef, basisobj, fdnames)
  return(stdfd)
}
