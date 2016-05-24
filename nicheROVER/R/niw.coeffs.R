niw.coeffs <-
function(X, lambda, kappa, Psi, nu) {
  d <- ncol(X)
  N <- nrow(X)
  if(missing(kappa)) kappa <- 0
  if(missing(Psi)) Psi <- 0
  if(missing(nu)) nu <- ncol(X)+1
  # sufficient statistics
  Xbar <- colMeans(X)
  S <- t(X)-Xbar
  S <- S %*% t(S)
  # posterior parameter values
  Psi2 <- Psi + S
  lambda2 <- N*Xbar
  if(kappa != 0) {
    Psi2 <- Psi2 + (N*kappa)/(N+kappa) * (Xbar-lambda) %*% t(Xbar-lambda)
    lambda2 <- lambda2 + kappa*lambda
  }
  lambda2 <- lambda2/(N+kappa)
  nu2 <- N+nu-(kappa==0)
  kappa2 <- N+kappa
  list(lambda = lambda2, kappa = kappa2, Psi = Psi2, nu = nu2)
}
