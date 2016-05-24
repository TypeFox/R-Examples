### Density function for the skew-Laplace distribution
dskewlap <- function(x, Theta, logPars = FALSE){
  if(length(Theta) != 3 )
    stop("parameter vector must contain 3 values")
  alpha <- if(logPars == TRUE) exp(Theta[1]) else Theta[1] 
  beta <- if(logPars == TRUE) exp(Theta[2]) else Theta[2]
  mu <- Theta[3] 
  ab <- alpha + beta
  belowMu <- (1/ab)*exp((x - mu)/alpha)
  aboveMu <- (1/ab)*exp((mu - x)/beta)
  skewlapDens <- ifelse(x <= mu, belowMu, aboveMu) 
  skewlapDens
} ## End of dskewlap()

### Distribution function for the skew-Laplace distribution
pskewlap <- function(q, Theta){
  if (length(Theta) != 3)
    stop("parameter vector must contain 3 values")
  alpha <- Theta[1] 
  beta <- Theta[2]
  mu <- Theta[3]
  ab <- alpha + beta

  if(alpha <= 0) stop("alpha must be positive")
  if(beta <= 0) stop("beta must be positive")

  belowMu <- (alpha/ab)*exp((q - mu)/alpha)
  aboveMu <- 1 - (beta/ab)*exp((mu - q)/beta)

  distFn <- ifelse(q < mu, belowMu, aboveMu)
  distFn
} ## End of pskewlap()


### Quantile function for the skew-Laplace distribution
qskewlap <- function(p, Theta){
  if (length(Theta) != 3)
    stop("parameter vector must contain 3 values")
  alpha <- Theta[1] 
  beta <- Theta[2]
  mu <- Theta[3]
  ab <- alpha + beta

  if(alpha <= 0) stop("alpha must be positive")
  if(beta <= 0) stop("beta must be positive")
  if(any(p < 0)||any(p > 1)) stop("p must be between 0 and 1")

  belowMu <- alpha*log(p*ab/alpha) + mu
  aboveMu <- mu - beta*log(ab*(1 - p)/beta)

  qFn <- ifelse(p < alpha/ab, belowMu, aboveMu)
  qFn
} ## End of qskewlap()

### Function to generate random observations from the
### skew Laplace distribution
rskewlap <- function(n, Theta){
  if (length(Theta) != 3)
    stop("parameter vector must contain 3 values")
  alpha <- Theta[1] 
  beta <- Theta[2]
  mu <- Theta[3]
  ab <- alpha + beta

  if(alpha <= 0) stop("alpha must be positive")
  if(beta <= 0) stop("beta must be positive")

  y <- rexp(n,1)
  probs <- c(alpha,beta)/ab
  signs <- sample(c(-1,1), n, replace = TRUE, prob = probs)
  mult <- ifelse(signs < 0, signs*alpha, signs*beta)
  x <- mult*y + mu
  x
} ## End of rskewlap()
