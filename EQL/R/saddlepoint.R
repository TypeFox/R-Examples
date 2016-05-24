.getFactor <- function(n, cum) {
  K      <- cum$K
  mu.inv <- cum$mu.inv
  kappa2 <- cum$kappa2
  domain <- cum$domain
  dens <- function(x) sqrt(n / kappa2(mu.inv(x))) *
    exp(n * (K(mu.inv(x)) - mu.inv(x) * x))
  return(1 / integrate(dens, domain[1], domain[2])$value)
}

saddlepoint <- function(x, n, cumulants, correct=TRUE, normalize=FALSE) {
  if (!check(cumulants)) {
    stop("Object 'cumulants' is not a properly defined 'cumulants' object!")
  }
  if (correct && cumulants$missing) {
    warning("The correction term cannot be used, for the higher ",
            "cumulants are not given!")
    correct <- FALSE
  }
  if (correct && normalize) {
    warning("The renormalized version does not use the correction term. ",
            "It will be skipped!")
    correct <- FALSE
  }
  z <- cumulants$mu.inv(x)
  kappa2 <- cumulants$kappa2(z)
  K <- cumulants$K(z)
  corTerm <- 1
  if (correct) {
    r3 <- cumulants$rho3(z)
    r4 <- cumulants$rho4(z)
    corTerm <- corTerm + (3 * r4 - 5 * r3 ^ 2) / (24 * n)
  }
  value <- sqrt(n / (2 * pi * kappa2)) * exp(n * K - n * z * x) * corTerm
  if (normalize) {
    value <- .getFactor(n, cumulants) * value
  }
  return(approximation(y=x, approx=value, n=n, type="mean",
                       approx.type="Saddlepoint"))
}
