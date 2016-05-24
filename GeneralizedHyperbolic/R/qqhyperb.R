### Q-Q plot for hyperbolic distribution
qqhyperb <- function(y, mu = 0, delta = 1, alpha = 1, beta = 0,
                     param = c(mu, delta, alpha, beta),
                     main = "Hyperbolic Q-Q Plot",
                     xlab = "Theoretical Quantiles",
                     ylab = "Sample Quantiles",
                     plot.it = TRUE, line = TRUE, ...) {

  qqghyp(y, param = c(param, 1), main = main, xlab = xlab,
         ylab = ylab, plot.it = plot.it, line = line, ...)
} ## End of qqhyperb()

### P-P plot for hyperbolic distribution
pphyperb <- function(y, mu = 0, delta = 1, alpha = 1, beta = 0,
                     param = c(mu, delta, alpha, beta),
                     main = "Hyperbolic P-P Plot",
                     xlab = "Uniform Quantiles",
                     ylab = "Probability-integral-transformed Data",
                     plot.it = TRUE, line = TRUE, ...) {

  ppghyp(y, param = c(param, 1), main = main, xlab = xlab,
         ylab = ylab, plot.it = plot.it, line = line, ...)
} ## End of pphyperb()


  
