EV.lnre.fzm <- function (obj, N=NA, ...)
{
  if (! inherits(obj, "lnre.fzm")) stop("argument must be object of class 'lnre.fzm'")
  if (missing(N)) stop("argument 'N' is required for 'lnre.fzm' objects")
  if (!(is.numeric(N) && all(N >= 0))) stop("argument 'N' must be non-negative integer")

  alpha <- obj$param$alpha
  A <- obj$param$A
  B <- obj$param$B
  C <- obj$param2$C

  if (obj$exact) {
    term1 <- (N ^ alpha) *
      (Igamma(1 - alpha, N * A, lower=FALSE) - Igamma(1 - alpha, N * B, lower=FALSE))
    term2 <- (1 - exp(- N * A)) / A^alpha - (1 - exp(- N * B)) / B^alpha
    (C / alpha) * (term1 + term2)
  }
  else {
    term1 <- (N ^ alpha) * Igamma(1 - alpha, N * A, lower=FALSE)
    term2 <- (A ^ (-alpha)) * (1 - exp(-N * A))
    (C / alpha) * (term1 + term2)
  }
}
