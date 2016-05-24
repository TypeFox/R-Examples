### Q-Q plot for normal inverse Gaussian distribution
qqnig <- function(y, mu = 0, delta = 1, alpha = 1, beta = 0,
                  param = c(mu, delta, alpha, beta),
                  main = "Normal inverse Gaussian Q-Q Plot",
                  xlab = "Theoretical Quantiles",
                  ylab = "Sample Quantiles",
                  plot.it = TRUE, line = TRUE, ...) {

  qqghyp(y, param = c(param, -1/2), main = main, xlab = xlab,
        ylab = ylab, plot.it = plot.it, line = line, ...)
} ## End of qqnig()

### P-P plot for normal inverse Gaussian distribution
ppnig <- function(y, mu = 0, delta = 1, alpha = 1, beta = 0,
                  param = c(mu, delta, alpha, beta),
                  main = "Normal inverse Gaussian P-P Plot",
                  xlab = "Uniform Quantiles",
                  ylab = "Probability-integral-transformed Data",
                  plot.it = TRUE, line = TRUE, ...) {

  ppghyp(y, param = c(param, -1/2), main = main, xlab = xlab,
        ylab = ylab, plot.it = plot.it, line = line, ...)
} ## End of ppnig()



