EVm.lnre.zm <- function (obj, m, N=NA, ...)
{
  if (! inherits(obj, "lnre.zm")) stop("argument must be object of class 'lnre.zm'")
  if (missing(N)) stop("argument 'N' is required for 'lnre.zm' objects")
  if (!(is.numeric(N) && all(N >= 0))) stop("argument 'N' must be non-negative integer")
  if (!(is.numeric(m) && all(m >= 1))) stop("argument 'm' must be positive integer")

  alpha <- obj$param$alpha
  B <- obj$param$B
  C <- obj$param2$C

  if (obj$exact) {
    factor <- exp(Igamma(m - alpha, N * B, log=TRUE) - Cgamma(m + 1, log=TRUE)) # = Igamma(m - alpha, N*B) / m!
    C * N^alpha * factor
  }
  else {
    factor <- exp(Cgamma(m - alpha, log=TRUE) - Cgamma(m + 1, log=TRUE)) # = Cgamma(m - alpha) / m!
    C * N^alpha * factor
  }
}
