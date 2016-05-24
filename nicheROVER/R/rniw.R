rniw <-
function(n, lambda, kappa, Psi, nu) {
  d <- length(lambda)
  Sigma <- rwish(n, Psi, nu, inv = TRUE)
  mu <- matrix(NA, n, d)
  colnames(mu) <- names(lambda)
  for(ii in 1:n) {
    mu[ii,] <- rmvnorm(1, mean = lambda, sigma = Sigma[,,ii]/kappa)
  }
  list(mu = mu, Sigma = Sigma)
}
