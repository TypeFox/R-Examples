#' Percentiles of ES distribution function for normally distributed geometric 
#' returns
#' 
#' Estimates the percentiles of ES distribution for normally distributed 
#' geometric returns, for specified confidence level and holding period using 
#' the theory of order statistics.
#' 
#' @param ... The input arguments contain either return data or else mean and 
#' standard deviation data. Accordingly, number of input arguments is either 5 
#' or 7. In case there 5 input arguments, the mean, standard deviation and 
#' number of samples is computed from return data. See examples for details.
#' 
#'  returns Vector of daily geometric return data
#' 
#'  mu Mean of daily geometric return data
#'  
#'  sigma Standard deviation of daily geometric return data
#'  
#'  n Sample size
#'  
#'  investment Size of investment
#'  
#'  perc Desired percentile
#'  
#'  cl ES confidence level and must be a scalar
#'  
#'  hp ES holding period and must be a a scalar
#'  
#' @return Percentiles of ES distribution function
#'  
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#'
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Estimates Percentiles of ES distribution
#'    data <- runif(5, min = 0, max = .2)
#'    LogNormalESDFPerc(returns = data, investment = 5, perc = .7, cl = .95, hp = 60)
#'    
#'    # Estimates Percentiles given mean, standard deviation and number of sambles of return data
#'    LogNormalESDFPerc(mu = .012, sigma = .03, n= 10, investment = 5, perc = .8, cl = .99, hp = 40)
#'
#'
#' @export
LogNormalESDFPerc <- function(...){
  # Determine if there are five or seven arguments, and ensure that arguments are read as intended
  if (nargs() < 5) {
    stop("Too few arguments")
  }
  if (nargs() == 6) {
    stop("Incorrect number of arguments")
  }
  if (nargs() > 7) {
    stop("Too many arguments")
  }
  args <- list(...)
  if (nargs() == 7) {
    mu <- args$mu
    investment <- args$investment
    cl <- args$cl
    perc <- args$sigma
    n <- args$n
    sigma <- args$sigma
    hp <- args$hp
  }
  if (nargs() == 5) {
    mu <- mean(args$returns)
    investment <- args$investment
    n <- max(dim(as.matrix(args$returns)))
    perc <- args$perc
    cl <- args$cl
    sigma <- sd(args$returns)
    hp <- args$hp
  }
  
  # Check that inputs have correct dimensions
  mu <- as.matrix(mu)
  mu.row <- dim(mu)[1]
  mu.col <- dim(mu)[2]
  if (max(mu.row, mu.col) > 1) {
    stop("Mean must be a scalar")
  }
  sigma <- as.matrix(sigma)
  sigma.row <- dim(sigma)[1]
  sigma.col <- dim(sigma)[2]
  if (max(sigma.row, sigma.col) > 1) {
    stop("Standard deviation must be a scalar")
  }
  n <- as.matrix(n)
  n.row <- dim(n)[1]
  n.col <- dim(n)[2]
  if (max(n.row, n.col) > 1) {
    stop("Number of observations in a sample must be a scalar")
  }
  if (n %% 1 != 0) {
    stop("Number of observations in a sample must be an integer.")
  }
  perc <- as.matrix(perc)
  perc.row <- dim(perc)[1]
  perc.col <- dim(perc)[2]
  if (max(perc.row, perc.col) > 1) {
    stop("Chosen percentile of the distribution must be a scalar")
  }
  cl <- as.matrix(cl)
  cl.row <- dim(cl)[1]
  cl.col <- dim(cl)[2]
  if (max(cl.row, cl.col) > 1) {
    stop("Confidence level must be a scalar")
  }
  hp <- as.matrix(hp)
  hp.row <- dim(hp)[1]
  hp.col <- dim(hp)[2]
  if (max(hp.row, hp.col) > 1) {
    stop("Holding period must be a scalar")
  }
  
  # Check that inputs obey sign and value restrictions
  if (sigma < 0) {
    stop("Standard deviation must be non-negative")
  }
  if (n < 0) {
    stop("Number of observations must be non-negative")
  }
  if (perc > 1){
    stop("Chosen percentile must not exceed 1")
  }
  if (perc <= 0){
    stop("Chosen percentile must be positive")
  }
  if (cl >= 1){
    stop("Confidence level must be less than 1")
  }
  if (cl <= 0){
    stop("Confidence level must be greater than 0")
  }
  if (hp <= 0){
    stop("Honding period must be greater than 0")
  }
  
  # Derive order statistic and ensure it is an integer
  w <- n * cl # Derive r-th order statistic
  r <- round(w) # Round r to nearest integer
  # Bisection routine
  a <- 0
  fa <- -Inf
  b <- 1
  fb <- Inf
  eps <- .Machine$double.eps
  while (b - a > eps * b) {
    x <- (a + b) / 2
    fx <- 1 - pbinom(r - 1, n, x) - perc
    if (sign(fx) == sign(fa)){
      a = x
      fa = fx
    } else {
      b = x
      fb = fx
    }
  }
  
  # VaR estimation
  VaR <- investment - exp( sigma[1,1] * sqrt(hp) %*% qnorm(1 - cl, 0, 1)  + mu[1,1] * hp %*% matrix(1,cl.row,cl.col) + log(investment)) # VaR
  
  # ES estimation
  cl <- x
  n <- 1000 # Number of slices into which tail is divided
  cl0 <- cl # Initial confidence level
  delta.cl <- (1 - cl) / n # Increment to confidence level as each slice is taken
  term <- VaR
  for (i in 1:(n-1)) {
    cl <- cl0 + i * delta.cl # Revised cl
    term <- term + investment - exp(sigma[1,1] * sqrt(hp) %*% qnorm(1 - cl, 0, 1)  + mu[1,1] * hp %*% matrix(1,cl.row,cl.col) + log(investment))
  }
  y <- term/n
  return(y)
  
}