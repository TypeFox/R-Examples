### Functions for the hyperbolic distribution
### Density of the hyperbolic distribution
dhyperb <- function(x, mu = 0, delta = 1, alpha = 1, beta = 0,
                    param = c(mu, delta, alpha, beta)) {

  if (length(param) != 4)
    stop("param vector must contain 4 values")

  param <- as.numeric(param)

  dghyp(x, param = c(param, 1))
} ## End of dhyperb()


### Cumulative distribution function of the hyperbolic distribution
### Now just calls pghyp
### DJS 28/7/10
### DJS 05/09/06
phyperb <- function(q, mu = 0, delta = 1, alpha = 1, beta = 0,
                    param = c(mu, delta, alpha, beta),
                    lower.tail = TRUE, subdivisions = 100,
                    intTol = .Machine$double.eps^0.25,
                    valueOnly = TRUE, ...){

  if (length(param) != 4)
    stop("param vector must contain 4 values")

  param <- as.numeric(param)

  pghyp(q, param = c(param, 1), lower.tail = lower.tail,
        subdivisions = subdivisions, intTol = intTol,
        valueOnly = valueOnly, ...)
} ## End of phyperb()

### qhyperb using breaks as for phyperb and splines as in original qhyperb
### Now just calls qghyp
### DJS 28/7/10
### DJS 06/09/06
qhyperb <- function(p, mu = 0, delta = 1, alpha = 1, beta = 0,
                    param = c(mu, delta, alpha, beta),
                    lower.tail = TRUE, method = c("spline", "integrate"),
                    nInterpol = 501, uniTol = .Machine$double.eps^0.25,
                    subdivisions = 100, intTol = uniTol, ...) {

  if (length(param) != 4)
   stop("param vector must contain 4 values")

  param <- as.numeric(param)

  qghyp(p, param = c(param,1), lower.tail = lower.tail, method = method,
        nInterpol = nInterpol, uniTol = uniTol,
        subdivisions = subdivisions, intTol = intTol, ...)
} # End of qhyperb()

### Function to generate random observations from a
### hyperbolic distribution using the
### mixing property of the generalized inverse
### Gaussian distribution and Dagpunar's algorithm
### for the generalized inverse Gaussian
rhyperb <- function(n, mu = 0, delta = 1, alpha = 1, beta = 0,
                    param = c(mu, delta, alpha, beta)) {

  if (length(param) != 4)
    stop("param vector must contain 4 values")

  param <- as.numeric(param)

  rghyp(n, param = c(param,1))
} ## End of rhyperb()

### Derivative of the density
ddhyperb <- function(x, mu = 0, delta = 1, alpha = 1, beta = 0,
                     param = c(mu, delta, alpha, beta)) {

  if (length(param) != 4)
    stop("param vector must contain 4 values")

  param <- as.numeric(param)

  ddghyp(x, param = c(param, 1))
} ## End of ddhyperb()

