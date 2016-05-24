niw.mom <-
function(lambda, kappa, Psi, nu) {
  d <- length(lambda)
  b <- nu-d
  mu.mean <- lambda
  Sigma.mean <- Psi/(b-1)
  mu.var <- Sigma.mean/kappa
  Sigma.var <- Psi %o% Psi
  Sigma.var <- 2 * Sigma.var + (b-1) * (aperm(Sigma.var, c(1,3,2,4)) +
                                        aperm(Sigma.var, c(1,4,3,2)))
  Sigma.var <- Sigma.var/(b*(b-1)^2*(b-3))
  list(mu = list(mean = mu.mean, var = mu.var),
       Sigma = list(mean = Sigma.mean, var = Sigma.var))
}
