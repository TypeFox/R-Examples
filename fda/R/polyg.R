polyg <- function(x, argvals, nderiv=0)
{
#  Evaluates the basis for a linear interpolant or its first derivative.
#  It calls function spline.des.
#  Arguments are as follows:
#  X      ... array of values at which the spline functions are to
#             evaluated
#  ARGVAL ... a STRICTLY INCREASING sequence of argument values.
#  NDERIV ... Either 0 or 1.  0 means only function values
#             are returned; 1 means derivative values are returned
#  Return is a matrix with length(X) rows and number of columns equal to
#             number of argument values

#  last modified 8 June 1999

  x <- as.vector(x)
  n <- length(x)

  if (!is.array(argvals)) argvals <- as.array(argvals)
  if (length(dim(argvals)) != 1) stop(
     'ARGVALS is not a vector or 1-dim. array.')
  if ( (max(x) > max(argvals)) || (min(x) < min(argvals)) ) stop(
     'ARGVALS do not span the values of X.')

  nargvals <- length(argvals)
  if (min(diff(argvals)) <= 0 ) stop(
     'Break-points are not strictly increasing')

  if (!(nderiv == 0 | nderiv == 1)) stop(
     'NDERIV is neither 0 nor 1.')
  derivs    <- rep(nderiv,n)
  nbasis <- length(argvals)

  knots <- c(argvals[1], argvals, argvals[nbasis])
  basismat <- splines::spline.des(knots, x, 2, derivs)$design

  return (basismat)
}
