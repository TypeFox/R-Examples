vgc.interp <- function (obj, N, m.max=0, allow.extrapolation=FALSE)
{
  if (! inherits(obj, "spc")) stop("first argument must be object of class 'spc'")
  if (attr(obj, "m.max") > 0) stop("cannot interpolate from incomplete frequency spectrum")
  if (attr(obj, "expected")) stop("cannot interpolate from expected frequency spectrum")
  if (! (is.numeric(N) && all(N >= 0))) stop("'N' must be a vector of non-negative integers")
  if (any(diff(N) < 0)) {
    warning("'N' must be increasing (data points have been reordered!)")
    N <- sort(N)
  }
  if (!missing(m.max) && !(length(m.max) == 1 && is.numeric(m.max) && 0 <= m.max && m.max <= 9))
    stop("'m.max' must be a single integer in the range 1 ... 9")

  N0 <- N(obj)
  if (any(N > N0) && !allow.extrapolation)
    stop("binomial extrapolation to N=", max(N), " from N0=", N0, " not allowed!")

  E.V <- EV(obj, N, allow.extrapolation=allow.extrapolation)
  if (m.max > 0) {
    E.Vm <- lapply(1:m.max,             # make list of growth curves for E[V_m(N)]
                   function (.m) EVm(obj, .m, N, allow.extrapolation=allow.extrapolation))
    vgc(N=N, V=E.V, Vm=E.Vm, expected=TRUE)
  }
  else {
    vgc(N=N, V=E.V, expected=TRUE) 
  }
}
