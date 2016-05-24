EVm.spc <- function (obj, m, N, allow.extrapolation=FALSE, ...)
{
  if (! inherits(obj, "spc")) stop("first argument must be object of class 'spc'")
  if (attr(obj, "m.max") > 0) stop("cannot interpolate from incomplete frequency spectrum")
  if (attr(obj, "expected")) stop("cannot interpolate from expected frequency spectrum")
  if (! (is.numeric(N) && all(N >= 0))) stop("'N' must be vector of non-negative numbers")
  if (! (is.numeric(m) && all(m >= 1))) stop("'m' must be vector of integers >= 1")
  if (!is.integer(m)) {
    if (any(m != floor(m))) stop("only integer values are allowed for 'm'")
    m <- as.integer(m)
  }
  m.length <- length(m)
  N.length <- length(N)
  n.items <- max(m.length, N.length)
  if (m.length == 1) {
    m <- rep(m, n.items)
  }
  else if (N.length == 1) {
    N <- rep(N, n.items)
  }
  else {                                # m.length > 1 && N.length > 1
    stop("either 'N' or 'm' must be a single value")
  }
  
  ## Baayen (2001), p. 65, Eq. (2.41)
  N0 <- N(obj)
  if (any(N > N0) && !allow.extrapolation)
    stop("binomial extrapolation to N=", max(N), " from N0=", N0, " not allowed!")
  E.Vm <- numeric(n.items)
  for (.i in 1:n.items) {
    .N <- N[.i]
    .m <- m[.i]
    .idx <- obj$m >= .m
    if (.N <= N0) {              # interpolation -> use binomial probabilities
      E.Vm[.i] <- sum( obj$Vm[.idx] * dbinom(.m, obj$m[.idx], .N / N0) )
    }
    else {                       # extrapolation -> compute alternating sum in the dumb way
      .f1 <- .N / N0
      .f2 <- 1 - .f1
      .k <- obj$m[.idx] # this is the summation variable k in Baayen's formula
      E.Vm[.i] <- sum( obj$Vm[.idx] * choose(.k, .m) * (.f1 ^ .m) * (.f2 ^ (.k-.m)) )
    }
  }

  E.Vm
}
