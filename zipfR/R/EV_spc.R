EV.spc <- function (obj, N, allow.extrapolation=FALSE, ...)
{
  if (! inherits(obj, "spc")) stop("first argument must be object of class 'spc'")
  if (attr(obj, "m.max") > 0) stop("cannot interpolate from incomplete frequency spectrum")
  if (attr(obj, "expected")) stop("cannot interpolate from expected frequency spectrum")
  if (! (is.numeric(N) && all(N >= 0))) stop("'N' must be vector of non-negative numbers")

  ## Baayen (2001), p. 65, Eq. (2.43)
  N0 <- N(obj)
  V0 <- V(obj)
  if (any(N > N0) && !allow.extrapolation)
    stop("binomial extrapolation to N=", max(N), " from N0=", N0, " not allowed!")
  n.samples <- length(N)
  E.V <- numeric(n.samples)
  for (.i in 1:n.samples) {
    .N <- N[.i]
    factor <- 1 - .N/N0
    E.V[.i] <- V0 - sum( obj$Vm * (factor ^ obj$m) )
  }

  E.V
}
