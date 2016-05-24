#' Plots t VaR against confidence level and holding period
#' 
#' Plots the VaR of a portfolio against confidence level and holding period 
#' assuming that P/L are t distributed, for specified confidence level and 
#'  holding period.
#' 
#' @param ... The input arguments contain either return data or else mean and 
#'  standard deviation data. Accordingly, number of input arguments is either 4 
#'  or 5. In case there 4 input arguments, the mean and standard deviation of 
#'  data is computed from return data. See examples for details.
#'  
#'  returns Vector of daily geometric return data
#' 
#'  mu Mean of daily geometric return data
#' 
#'  sigma Standard deviation of daily geometric return data
#' 
#'  df Number of degrees of freedom in the t distribution
#' 
#'  cl VaR confidence level and must be a vector
#' 
#'  hp VaR holding period and must be a vector
#'  
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#'
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Plots VaR against confidene level given geometric return data
#'    data <- runif(5, min = 0, max = .2)
#'    tVaRPlot3D(returns = data, df = 6, cl = seq(.85,.99,.01), hp = 60:90)
#'    
#'    # Computes VaR against confidence level given mean and standard deviation of return data
#'    tVaRPlot3D(mu = .012, sigma = .03, df = 6, cl = seq(.85,.99,.02), hp = 40:80)
#'
#'
#' @export
tVaRPlot3D <- function(...){
  if (nargs() < 4) {
    stop("Too few arguments")
  }
  if (nargs() > 5) {
    stop("Too many arguments")
  }
  args <- list(...)
  if (nargs() == 5) {
    mu <- args$mu
    df <- args$df
    cl <- args$cl
    sigma <- args$sigma
    hp <- args$hp
  }
  if (nargs() == 4) {
    mu <- mean(args$returns)
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
    stop("Confidence level must be a vector")
  }
  hp <- as.matrix(hp)
  hp.row <- dim(hp)[1]
  hp.col <- dim(hp)[2]
  if (min(hp.row, hp.col) > 1) {
    stop("Holding period must be a vector")
  }
  df <- as.matrix(df)
  df.row <- dim(df)[1]
  df.col <- dim(df)[2]
  if (max(df.row, df.col) > 1) {
    stop("Number of degrees of freedom must be a scalar")
  }
  
  # Check that cl is read as row vector
  if (cl.row > cl.col) {
    cl <- t(cl)
  }
  # Check that hp is read as column vector
  if (hp.row > hp.col) {
    hp <- t(hp)
  }
  
  # Check that inputs obey sign and value restrictions
  if (sigma < 0) {
    stop("Standard deviation must be non-negative")
  }
  if (df < 3) {
    stop("Number of degrees of freedom must be at least 3 for first two moments of distribution to be defined")
  }
  if (max(cl) >= 1){
    stop("Confidence level(s) must be less than 1")
  }
  if (min(cl) <= 0){
    stop("Confidence level(s) must be greater than 0")
  }
  if (min(hp) <= 0){
    stop("Holding period(s) must be greater than 0")
  }
  
  # VaR estimation
  cl.row <- dim(cl)[1]
  cl.col <- dim(cl)[2]
  # VaR estimation
  VaR <- (-sigma[1,1] * sqrt(t(hp)) %*% sqrt((df - 2) / df) %*% qt(1 - cl, df)) + (- mu[1,1] * t(hp) %*% matrix(1, cl.row, cl.col)) # VaR
  # Plotting
  persp(x=cl, y=hp, t(VaR), xlab = "Confidence Level", 
        ylab = "Holding Period", zlab = "VaR", 
        main = "t VaR against Confidence Level and Holding Period")
  
}
