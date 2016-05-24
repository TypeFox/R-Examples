### Density function for the skew-Laplace distribution
dskewlap <- function(x, mu = 0, alpha = 1, beta = 1,
                     param = c(mu, alpha, beta),
                     logPars = FALSE) {

  if (length(param) != 3)
    stop("param vector must contain 3 values")

  mu <- param[1]
  alpha <- if(logPars == TRUE) exp(param[2]) else param[2]
  beta <- if(logPars == TRUE) exp(param[3]) else param[3]
  ab <- alpha + beta

  belowMu <- (1 / ab) * exp((x - mu) / alpha)
  aboveMu <- (1 / ab) * exp((mu - x) / beta)
  skewlapDens <- ifelse(x <= mu, belowMu, aboveMu)

  skewlapDens
} ## End of dskewlap()


### Distribution function for the skew-Laplace distribution
pskewlap <- function(q, mu = 0, alpha = 1, beta = 1,
                     param = c(mu, alpha, beta)) {

  if (length(param) != 3)
    stop("param vector must contain 3 values")

  mu <- param[1]
  alpha <- param[2]
  beta <- param[3]
  ab <- alpha + beta

  if (alpha <= 0)
    stop("alpha must be positive")

  if (beta <= 0)
    stop("beta must be positive")

  belowMu <- (alpha / ab) * exp((q - mu) / alpha)
  aboveMu <- 1 - (beta / ab) * exp((mu - q) / beta)

  distFn <- ifelse(q < mu, belowMu, aboveMu)
  distFn
} ## End of pskewlap()


### Quantile function for the skew-Laplace distribution
qskewlap <- function(p, mu = 0, alpha = 1, beta = 1,
                     param = c(mu, alpha, beta)) {

  if (length(param) != 3)
    stop("param vector must contain 3 values")

  mu <- param[1]
  alpha <- param[2]
  beta <- param[3]
  ab <- alpha + beta

  if (alpha <= 0)
    stop("alpha must be positive")

  if (beta <= 0)
    stop("beta must be positive")

  if (any(p < 0) | any(p > 1))
    stop("p must be between 0 and 1")

  belowMu <- alpha * log(p * ab / alpha) + mu
  aboveMu <- mu - beta * log(ab * (1 - p) / beta)

  qFn <- ifelse(p < alpha / ab, belowMu, aboveMu)
  qFn
} ## End of qskewlap()


### Function to generate random observations from the
### skew Laplace distribution
rskewlap <- function(n, mu = 0, alpha = 1, beta = 1,
                     param = c(mu, alpha, beta)) {

  if (length(param) != 3)
    stop("param vector must contain 3 values")

  mu <- param[1]
  alpha <- param[2]
  beta <- param[3]
  ab <- alpha + beta

  if (alpha <= 0)
    stop("alpha must be positive")

  if (beta <= 0)
    stop("beta must be positive")

  y <- rexp(n, 1)
  probs <- c(alpha, beta) / ab
  signs <- sample(c(-1, 1), n, replace = TRUE, prob = probs)
  mult <- ifelse(signs < 0, signs * alpha, signs * beta)
  x <- mult * y + mu
  x
} ## End of rskewlap()
