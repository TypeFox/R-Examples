#' Percentiles of VaR distribution function for normally distributed P/L
#' 
#' Estimates the percentile of VaR distribution function for normally distributed P/L, using the theory of order statistics.
#' 
#' @param ... The input arguments contain either return data or else mean and 
#' standard deviation data. Accordingly, number of input arguments is either 4 
#' or 6. In case there 4 input arguments, the mean, standard deviation and number of observations of 
#' data are computed from returns data. See examples for details.
#' 
#' returns Vector of daily geometric return data
#' 
#'  mu Mean of daily geometric return data

#'  sigma Standard deviation of daily geometric return data
#' 
#'  n Sample size
#' 
#'  perc Desired percentile
#' 
#'  cl VaR confidence level and must be a scalar
#' 
#'  hp VaR holding period and must be a a scalar
#'  
#' @return Percentiles of VaR distribution function and is scalar
#'  
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#'
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Estimates Percentiles of VaR distribution
#'    data <- runif(5, min = 0, max = .2)
#'    NormalVaRDFPerc(returns = data, perc = .7, cl = .95, hp = 60)
#'    
#'    # Estimates Percentiles of VaR distribution
#'    NormalVaRDFPerc(mu = .012, sigma = .03, n= 10, perc = .8, cl = .99, hp = 40)
#'
#'
#' @export
NormalVaRDFPerc <- function(...){
  # Determine if there are four or six arguments, and ensure that arguments are read as intended
  if (nargs() < 4) {
    stop("Too few arguments")
  }
  if (nargs() == 5) {
    stop("Incorrect number of arguments")
  }
  if (nargs() > 6) {
    stop("Too many arguments")
  }
  args <- list(...)
  if (nargs() == 6) {
    mu <- args$mu
    cl <- args$cl
    perc <- args$sigma
    n <- args$n
    sigma <- args$sigma
    hp <- args$hp
  }
  if (nargs() == 4) {
    mu <- mean(args$returns)
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
    stop("Number of observations in a sample must be an integer")
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
    stop("Confidence level(s) must be less than 1")
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
  y <- qnorm(x, -mu * hp, sigma*sqrt(hp)) # Value of percentile; not normal VaR formula
  
  return(y)
}
