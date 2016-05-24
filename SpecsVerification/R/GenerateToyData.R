GenerateToyData <- function(N = 20, mu.y = 0, s.s = 7, s.eps = 6, mu.x = 0, beta = 0.2, s.eta = 8, K = 10, mu.x.ref = NA, beta.ref = NA, s.eta.ref = NA, K.ref = NA) {
  s <- rnorm(N, 0, s.s)
  y <- mu.y + s + rnorm(N, 0, s.eps)
  x1 <- mu.x + beta * s + matrix(rnorm(N*K, 0, s.eta), N, K)
  if (!any(is.na(c(mu.x.ref, beta.ref, s.eta.ref, K.ref)))) {
    x2 <- mu.x.ref + beta.ref * s + matrix(rnorm(N*K.ref, 0, s.eta.ref), N, K.ref)
    ret <- list(obs=y, ens=x1, ens.ref=x2)
  } else {
    ret <- list(obs=y, ens=x1)
  }
  return(ret)
}

