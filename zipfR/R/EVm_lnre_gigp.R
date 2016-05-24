EVm.lnre.gigp <- function (obj, m, N=NA, ...)
{
  if (! inherits(obj, "lnre.gigp")) stop("argument must be object of class 'lnre.gigp'")
  if (missing(N)) stop("argument 'N' is required for 'lnre.gigp' objects")
  if (!(is.numeric(N) && all(N >= 0))) stop("argument 'N' must be non-negative integer")
  if (!(is.numeric(m) && all(m >= 1))) stop("argument 'm' must be positive integer")

  gamma <- obj$param$gamma
  b <- obj$param$B                      # use original notation from Baayen (2001)
  c <- obj$param$C
  Z <- obj$param2$Z

  ## TODO: re-implement using recurrence relation for V_m / alpha_m (Baayen 2001, p. 91)
  ##       because factor2 becomes very small and factor3 very large
  ## (probably requires C code for good performance because of recursive nature of code)
  term <- 1 + N/Z
  factor1 <- 2 * Z / (b * besselK(b, gamma+1) * term^(gamma/2)) # Baayen (2001), p. 90
  factor2 <- ( b * N / (2 * Z * sqrt(term)) )^(m) / Cgamma(m+1)
  factor3 <- besselK(b * sqrt(term), m + gamma)
  factor1 * factor2 * factor3
}
