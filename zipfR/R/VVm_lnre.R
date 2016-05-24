VVm.lnre <- function (obj, m, N=NA, ...)
{
  if (! inherits(obj, "lnre")) stop("first argument must be object of class 'lnre'")
  if (missing(N)) stop("argument 'N' is required for 'lnre' objects")
  if (!(is.numeric(N) && all(N >= 0))) stop("argument 'N' must be non-negative number")
  if (!(is.numeric(m) && all(m >= 1))) stop("argument 'm' must be positive integer")
  
  factor <- exp(lchoose(2*m, m) - 2*m * log(2)) # = choose(2*m, m) / 2^(2*m)
  EVm(obj, m, N) - factor * EVm(obj, 2*m, 2*N)
}
