# Function from package 'fda' (c) 2014

fourier <- function(x, nbasis = n, period = span, nderiv = 0)
{
  #  Computes the NDERIV derivative of the Fourier series basis
  #    for NBASIS functions with period PERIOD, these being evaluated
  #    at values in vector X
  #  Returns an N by NBASIS matrix of function values
  #  Note:  The number of basis functions always odd.  If the argument
  #   NBASIS is even, it is increased by one.

  #  last modified 23 February 2007 by Spencer Graves
  #  Previously modified 15 December 2005

  #  check x and set up range

  xNames <- names(x)
#
  x      <- as.vector(x)
  n      <- length(x)
  onen   <- rep(1,n)
  xrange <- range(x)
  span   <- xrange[2] - xrange[1]

  #  check period and set up omega

  if (period <= 0) stop("PERIOD not positive.")
  omega  <- 2*pi/period
  omegax <- omega*x

  #  check nbasis

  if (nbasis <= 0) stop("NBASIS not positive")

  #  check nderiv

  if (nderiv <  0) stop("NDERIV is negative.")

  #  if nbasis is even, add one

  if ((nbasis %/% 2)*2 == nbasis) nbasis <- nbasis + 1

  #  compute basis matrix

  basismat <- matrix(0,n,nbasis)
  if (nderiv == 0) {
    #  The fourier series itself is required.
    basismat[,1] <- 1/sqrt(2)
    if(nbasis>1){
      j    <- seq(2,nbasis-1,2)
      k    <- j/2
      args <- outer(omegax,k)
      basismat[,j]   <- sin(args)
      basismat[,j+1] <- cos(args)
    }
  } else {
    #  The derivative of order nderiv is required.
    basismat[,1] <- 0.0
    if(nbasis>1){
      if (nderiv == (nderiv %/% 2)*2) {
        mval  <- nderiv/2
        ncase <- 1
      } else {
        mval <- (nderiv-1)/2
        ncase <- 2
      }
      j    <- seq(2,nbasis-1,2)
      k    <- j/2
      fac  <- outer(onen,((-1)^mval)*(k*omega)^nderiv)
      args <- outer(omegax,k)
      if (ncase == 1) {
        basismat[,j]   <-  fac * sin(args)
        basismat[,j+1] <-  fac * cos(args)
      } else {
        basismat[,j]   <-  fac * cos(args)
        basismat[,j+1] <- -fac * sin(args)
      }
    }
  }
  basismat <- basismat/sqrt(period/2)
#
  fNames <- "const"
  n2 <- (nbasis %/%2)
  if(n2>0){
    SC <- outer(c("sin", "cos"), 1:n2, paste, sep="")
    fNames <- c(fNames, as.vector(SC))
  }
#
  dimnames(basismat) <- list(xNames, fNames)
#
  return(basismat)
}
