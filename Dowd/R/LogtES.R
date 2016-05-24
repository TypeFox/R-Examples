#' ES for t distributed geometric returns
#' 
#' Estimates the ES of a portfolio assuming that geometric returns are 
#' Student-t distributed, for specified confidence level and holding period.
#' 
#' @param ... The input arguments contain either return data or else mean and 
#'  standard deviation data. Accordingly, number of input arguments is either 5 
#'  or 6. In case there 5 input arguments, the mean and standard deviation of 
#'  data is computed from return data. See examples for details.
#'  
#'  returns Vector of daily geometric return data
#' 
#'  mu Mean of daily geometric return data
#' 
#'  sigma Standard deviation of daily geometric return data
#' 
#'  investment Size of investment
#' 
#'  df Number of degrees of freedom in the t distribution
#' 
#'  cl VaR confidence level
#' 
#'  hp VaR holding period
#'  
#' @return Matrix of ES whose dimension depends on dimension of hp and cl. If 
#' cl and hp are both scalars, the matrix is 1 by 1. If cl is a vector and hp is
#'  a scalar, the matrix is row matrix, if cl is a scalar and hp is a vector, 
#'  the matrix is column matrix and if both cl and hp are vectors, the matrix 
#'  has dimension length of cl * length of hp.
#'  
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#'
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Computes ES given geometric return data
#'    data <- runif(5, min = 0, max = .2)
#'    LogtES(returns = data, investment = 5, df = 6, cl = .95, hp = 90)
#'    
#'    # Computes ES given mean and standard deviation of return data
#'    LogtES(mu = .012, sigma = .03, investment = 5, df = 6, cl = .95, hp = 90)
#'
#'
#' @export
LogtES <- function(...){
  if (nargs() < 5) {
    stop("Too few arguments")
  }
  if (nargs() > 6) {
    stop("Too many arguments")
  }
  args <- list(...)
  if (nargs() == 6) {
    mu <- args$mu
    investment <- args$investment
    df <- args$df
    cl <- args$cl
    sigma <- args$sigma
    hp <- args$hp
  }
  if (nargs() == 5) {
    mu <- mean(args$returns)
    investment <- args$investment
    df <- args$df
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
  cl <- as.matrix(cl)
  cl.row <- dim(cl)[1]
  cl.col <- dim(cl)[2]
  if (min(cl.row, cl.col) > 1) {
    stop("Confidence level must be a scalar or a vector")
  }
  hp <- as.matrix(hp)
  hp.row <- dim(hp)[1]
  hp.col <- dim(hp)[2]
  if (min(hp.row, hp.col) > 1) {
    stop("Holding period must be a scalar or a vector")
  }
  
  # Check that cl and hp are read as row and column vectors respectively
  if (cl.row > cl.col) {
    cl <- t(cl)
  }
  if (hp.row > hp.col) {
    hp <- t(hp)
  }
  
  # Check that inputs obey sign and value restrictions
  if (sigma < 0) {
    stop("Standard deviation must be non-negative")
  }
  if (max(cl) >= 1){
    stop("Confidence level(s) must be less than 1")
  }
  if (min(cl) <= 0){
    stop("Confidence level(s) must be greater than 0")
  }
  if (min(hp) <= 0){
    stop("Holding Period(s) must be greater than 0")
  }
  # VaR estimation
  cl.row <- dim(cl)[1]
  cl.col <- dim(cl)[2]
  VaR <- investment - exp(((df - 2) / df) * sigma %*% sqrt(t(hp)) %*% qt(1 - cl, df)
                        + mu * t(hp) %*% matrix(1, cl.row, cl.col) + log(investment)) # VaR
  # ES estimation
  n <- 1000 # Number of slices into which tail is divided
  cl0 <- cl # Initial confidence level
  delta.cl <- (1 - cl) / n # Increment to confidence level as each slice is taken
  term <- VaR
  for (i in 1:(n-1)) {
    cl <- cl0 + i * delta.cl # Revised cl
    term <- term + investment - exp(((df - 2) / df) * sigma %*% sqrt(t(hp)) %*% 
                                      qt(1 - cl, df) + mu * t(hp) %*% 
                                      matrix(1, cl.row, cl.col) + log(investment))
  }
  y <- term/n
  return (y)
}