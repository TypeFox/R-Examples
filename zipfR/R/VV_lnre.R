VV.lnre <- function (obj, N=NA, ...)
{
  if (! inherits(obj, "lnre")) stop("first argument must be object of class 'lnre'")
  if (missing(N)) stop("argument 'N' is required for 'lnre' objects")
  if (!(is.numeric(N) && N >= 0)) stop("argument 'N' must be non-negative integer")
  
  EV(obj, 2*N) - EV(obj, N)
}
