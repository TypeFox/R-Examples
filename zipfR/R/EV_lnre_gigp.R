EV.lnre.gigp <- function (obj, N=NA, ...)
{
  if (! inherits(obj, "lnre.gigp")) stop("argument must be object of class 'lnre.gigp'")
  if (missing(N)) stop("argument 'N' is required for 'lnre.gigp' objects")
  if (!(is.numeric(N) && all(N >= 0))) stop("argument 'N' must be non-negative integer")

  gamma <- obj$param$gamma
  b <- obj$param$B                      # use original notation from Baayen (2001)
  c <- obj$param$C
  Z <- obj$param2$Z

  term <- 1 + N/Z
  factor1 <- 2 * Z / b                  # Baayen (2001), p. 90
  factor2 <- besselK(b, gamma) / besselK(b, gamma+1) # factor1 * factor2 = S
  factor3 <- besselK(b * sqrt(term), gamma) / (term^(gamma/2) * besselK(b, gamma))
  factor1 * factor2 * (1 - factor3)
}
