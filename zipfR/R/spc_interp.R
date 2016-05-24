spc.interp <- function (obj, N, m.max=max(obj$m), allow.extrapolation=FALSE)
{
  if (! inherits(obj, "spc")) stop("first argument must be object of class 'spc'")
  if (attr(obj, "m.max") > 0) stop("cannot interpolate from incomplete frequency spectrum")
  if (attr(obj, "expected")) stop("cannot interpolate from expected frequency spectrum")
  if (! (is.numeric(N) && all(N >= 0) && length(N) == 1))
    stop("'N' must be a single non-negative integer")
  N0 <- N(obj)
  if (any(N > N0) && !allow.extrapolation)
    stop("binomial extrapolation to N=", max(N), " from N0=", N0, " not allowed!")
  if (m.max < 1) stop("'m.max' must be >= 1")

  E.Vm <- EVm(obj, 1:m.max, N, allow.extrapolation=allow.extrapolation)
  if (missing(m.max)) { # don't use EV() for full spectrum, may be inconsistent!
    spc(Vm=E.Vm, m=1:m.max, expected=TRUE)
  }
  else {
    E.V <- EV(obj, N, allow.extrapolation=allow.extrapolation)
    spc(Vm=E.Vm, m=1:m.max, V=E.V, N=N, m.max=m.max, expected=TRUE)
  }
}
