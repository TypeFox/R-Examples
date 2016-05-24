### Function to calculate an effective range for the distribution function
### or for the density function of the Hyperbolic Distribution
### DJS 8/09/06
hyperbCalcRange <- function(mu = 0, delta = 1, alpha = 1, beta = 0,
                            param = c(mu, delta, alpha, beta),
                            tol = 10^(-5), density = TRUE, ...) {

  if (length(param) != 4)
    stop("param vector must contain 4 values")

  param <- as.numeric(param)

  ghypCalcRange(param = c(param, 1), tol = tol, density = density, ...)
} ## End of hyperbCalcRange()
