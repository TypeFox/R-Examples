expon <- function (x, ratevec=1, nderiv=0)
{
#  This computes values of the exponentials, or their derivatives.
#  RATEVEC is a vector containing the rate constants, or mulipliers of X
#    in the exponent of e.
#  The default is the exponential function.
#  Arguments are as follows:
#  X       ... array of values at which the polynomials are to
#             evaluated
#  RATEVEC ... a vector containing the rate constants, or mulipliers of X
#              in the exponent of e.
#  NDERIV  ... order of derivative.  0 means only function values
#             are returned.
#  Return is a matrix with length(X) rows and NRATE columns containing
#  the values of the exponential functions or their derivatives.

#  last modified 5 December 2001

  x <- as.vector(x)
  n <- length(x)
  nrate <- length(as.vector(ratevec))
  expval <- matrix(0,n,nrate)
  for (irate in 1:nrate) {
    rate <- ratevec[irate]
    expval[,irate] <- rate^nderiv * exp(rate*x)
  }
  return (expval)

}
