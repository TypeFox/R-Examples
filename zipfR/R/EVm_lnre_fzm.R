EVm.lnre.fzm <- function (obj, m, N=NA, ...)
{
  if (! inherits(obj, "lnre.fzm")) stop("argument must be object of class 'lnre.fzm'")
  if (missing(N)) stop("argument 'N' is required for 'lnre.fzm' objects")
  if (!(is.numeric(N) && all(N >= 0))) stop("argument 'N' must be non-negative integer")
  if (!(is.numeric(m) && all(m >= 1))) stop("argument 'm' must be positive integer")

  alpha <- obj$param$alpha
  A <- obj$param$A
  B <- obj$param$B
  C <- obj$param2$C

  if (obj$exact) {
    term1 <- exp(Igamma(m - alpha, N * A, lower=FALSE, log=TRUE) - Cgamma(m + 1, log=TRUE))
    # term1 = Igamma(m - alpha, N * A, lower=FALSE) / m!
    term2 <- exp(Igamma(m - alpha, N * B, lower=FALSE, log=TRUE) - Cgamma(m + 1, log=TRUE))
    # term1 = Igamma(m - alpha, N * B, lower=FALSE) / m!
    C * N^alpha * (term1 - term2)
  }
  else {
    factor <- exp(Igamma(m - alpha, N * A, lower=FALSE, log=TRUE) - Cgamma(m + 1, log=TRUE))
    # factor = Igamma(m - alpha, N * A, lower=FALSE) / m!
    C * N^alpha * factor
  }
}
