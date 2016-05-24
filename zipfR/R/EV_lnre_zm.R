EV.lnre.zm <- function (obj, N=NA, ...)
{
  if (! inherits(obj, "lnre.zm")) stop("argument must be object of class 'lnre.zm'")
  if (missing(N)) stop("argument 'N' is required for 'lnre.zm' objects")
  if (!(is.numeric(N) && all(N >= 0))) stop("argument 'N' must be non-negative integer")

  alpha <- obj$param$alpha
  B <- obj$param$B
  C <- obj$param2$C

  if (obj$exact) {
    (C * N^alpha / alpha) *
      (Igamma(1 - alpha, N * B) - (1 - exp(- N * B)) * (N * B)^(-alpha))
  }
  else {
    C * N^alpha * Cgamma(1 - alpha) / alpha
  }
}
