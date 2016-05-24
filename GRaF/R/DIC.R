DIC <-
function(object) {
  p <- pred(object$x, object, std = T, same = T)
  n <- length(object$obsy)

  # define log likelihood function
  lik <- function(x, y) d0(x, y)
  # deviance on mean of latent variable
  dtheta <- -2 * sum(object$wt * lik(p[, 1], object$obsy))
  # upper and lower bounds for numerical integration
  up <- p[, 1] + 10 * p[, 2]
  low <- p[, 1] - 10 * p[, 2]
  # log likelihood function to integrate over
  llik <- function(x, y, mn, std, wt) wt * lik(x, y) * dnorm(x, mn, std)
  # function to apply integration to a datapoint
  foo <- function(i, y, mu, sd, low, up, wt) {
    integrate(llik, low[i], up[i], y[i], mu[i], sd[i], wt[i])$val
  }

  # loop over datapoints, sum log-likelihoods and turn into mean deviance
  dbar <- -2 * sum(sapply(1:n, foo, object$obsy, p[, 1], p[, 2], low, up, object$wt))
 
  pD <- dbar - dtheta
  DIC <- dbar + pD
  c(DIC = DIC, pD = pD)
}