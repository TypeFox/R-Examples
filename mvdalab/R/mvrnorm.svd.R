mvrnorm.svd <- function (n = 1, mu = NULL, Sigma = NULL, tol = 1e-06, empirical = FALSE, Dist = "normal", 
                          skew = 5, skew.mean = 0, skew.sd = 1, poisson.mean = 5) {
input.vars <- mvrnormBase.svd(n = n, rep(0, times = ncol(Sigma)), Sigma = Sigma, empirical = empirical)
pvars <- pnorm(input.vars)
if(Dist == "poisson") {
  output <- qpois(pvars, poisson.mean)
} else if(Dist == "exp") {
  output <- qexp(pvars)
} else if(Dist == "skewnorm") {
  output <- qsn(pvars, xi = skew.mean, omega = skew.sd, alpha = skew); output <- matrix(output, ncol = ncol(Sigma), nrow = n)
} else if(Dist == "normal") {
  output <- mvrnormBase.svd(n = n, mu = mu, Sigma = Sigma, empirical = empirical)
} else stop("No distributions chosen")
output
}