#' Estimates the spectral risk measure of a portfolio 
#' 
#' Function estimates the spectral risk measure of a portfolio 
#' assuming losses are normally distributed, assuming exponential weighting
#' function with specified gamma.
#' 
#' @param mu Mean losses
#' @param sigma Standard deviation of losses
#' @param gamma Gamma parameter in exponential risk aversion
#' @param number.of.slices Number of slices into which density function is divided
#'
#' @return Estimated spectral risk measure
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Generates 95% confidence intervals for normal VaR for given parameters
#'    NormalSpectralRiskMeasure(0, .5, .8, 20)
#'    
#' @export
NormalSpectralRiskMeasure <- function(mu, sigma, gamma, number.of.slices){
  # Check that inputs obey sign and value restrictions
  if (sigma < 0) {
    stop("Standard deviation must be n.on-negative")
  }
  if (min(gamma) <= 0) {
    stop("Gamma must be greater than 0")
  }
  
  n <- number.of.slices
  
  # Crude (weighted average quantile) evstimate of risk measure

  crude.estimate.of.risk.measure <- crude.estimate.of.spectral.risk.measure(mu, sigma, gamma, n)
  crude.halving.error <- crude.estimate.of.risk.measure - crude.estimate.of.spectral.risk.measure(mu, sigma, gamma, n/2)
  # Trapezoidal rule estimate of risk measure
  trapezoidal.estimate <- trapezoidal.quadrature.estimate(mu, sigma, gamma, n)
  trapezoidal.halving.error <- (1/3) * trapezoidal.estimate - trapezoidal.quadrature.estimate(mu, sigma, gamma, n/2)
  # Simpson's rule estimate of risk measure
  simpsons.estimate <- simpsons.quadrature.estimate(mu, sigma, gamma, n)
  simpsons.halving.error <- (1/15) * (simpsons.estimate-simpsons.quadrature.estimate(mu, sigma, gamma, n/2))
  print(paste("Crude Estimate Of Risk Measure:", crude.estimate.of.risk.measure))
  print(paste("Crude Halving Error:", crude.halving.error))
  print(paste("Trapezoidal Estimate:",trapezoidal.estimate ))
  print(paste("Trapezoidal Halving Error:", trapezoidal.halving.error))
  print(paste("Simpsons Estimate:", simpsons.estimate))
  print(paste("Simpsons Halving Error:", simpsons.halving.error))
  
}

crude.estimate.of.spectral.risk.measure <- function(mu, sigma, gamma, n) {
  # Applies crude average approach to estimate exponential spectral risk measure
  # Input arguments:
  #    mu : mean losses
  #    sigma : std. losses
  #    gamma : gamma weight in exponential spectral risk aversion function
  #    n : number of slices
  p <- seq(1/n, (n-1)/n, 1/n)
  product <- double(n-1)
  phi <- double(n-1)
  VaR <- double(n-1)
  for (i in 1:(n-1)) {
    VaR[i] <- mu + sigma * qnorm(p[i], 0, 1) # VaRs
    phi[i]=exp(-(1-p[i])/gamma)/(gamma*(1-exp(-1/gamma))); # Weights
    product[i]=VaR[i]*phi[i]; # Weighed VaR
  }
  y <- sum(product) / (n - 1) # Crude estimate of exponential spectral risk measure
  return(y)
}


trapezoidal.quadrature.estimate<- function(mu, sigma, gamma, n) {
  # Applies trapezioidal rule to estimateintegral of f( x) numerifcally using trapezoidal rule with given n, where f(x) is the fuction in te exponentnial spectrahl risk measure.
  # Input parameters:
  #    mu : mean losses
  #    sigma : standard losses
  #    gamma : gamma wight in exponential sepectral risk aversion function
  #    n : number of slics
  a <- 1/n
  b <- (n-1)/n # Limits of integration, bearing in mind we wish to avoid limits of 0 and 1 because inverses may not be not defined
  h <- (b-a)/(n-1)
  p <- double(n)
  for (i in 1:n) {
    p[i] <- a + (i - 1) * h
  }
  w <- double(n)
  w[1] <- h/2  # Initial trap weights
  w[n] <- h/2  # Other trap weights
  for (i in 2:(n-1)) {
    w[i] <- h
  }
  # Specify f(x)
  phi <- double(n)
  VaR <- double(n)
  f <- double(n)
  for (i in 1:n) {
    VaR[i] <- mu + sigma * qnorm(p[i], 0, 1) # VaRs
    phi[i] <- exp(-(1-p[i]) /gamma)/(gamma * (1-exp(-1/gamma))) # Spectral weights in  risk measure
    f[i] <- VaR[i] * phi[i] # f(i), weighted VaR
  }
  y <- t(as.matrix(w)) %*% as.matrix(f)
  return(y)
}

simpsons.quadrature.estimate <- function(mu, sigma, gamma, n) {
  # Function applies Simpson's rule to estimate intaegral of f(x) numerically
  # using Simpson's rule with given n, where f(x) is the function in the 
  # exponential spectral risk measure.
  # Input arguments:
  #    mu : mean losses
  #    sigma : std losses
  #    gamma : gamma weight in exponential spectral risk aversion function
  #    n : number of slices. NB: must be even
  n <- n - 1 # Convert to odd for purposes of algorithm
  a <- 1/n
  b <- (n - 1) / n # Limits of integration, bearing in mind we wish to avoid limits of 0 and 1 because inverses may not be not defined
  h <- (b - a) / (n - 1) # Increment
  p <- double(n)
  for (i in 1:n) { # Domain of integration, x
    p[i] <- a + (i - 1) * h
  }
  # Simpson's rule weights
  a[1] <- h/3
  w <- double(n-1)
  w[n] <- h/3 # Initial trap weights
  for (i in seq(2, (n-1), 2)) { # odd trap weights
    w[i] <- 4 * h / 3
  }
  # Specify f(x)
  VaR <- double(n)
  phi <- double(n)    
  f <- double(n)
  for (i in 1:n) {
    VaR[i] <- mu + sigma *qnorm(p[i], 0, 1) # VaRs
    phi[i] <- exp( - (1 - p[i]) / gamma) / (gamma * (1-exp(-1/gamma)))
    # Spectral weights in risk measure
    f[i] <- VaR[i] * phi[i] # f[i], weighted VaR
  }
  y <- t(as.matrix(w)) %*% as.matrix(f) 
  return(y)
}