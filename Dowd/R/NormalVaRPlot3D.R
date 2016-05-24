#' Plots normal VaR in 3D against confidence level and holding period
#' 
#' Plots the VaR of a portfolio against confidence level and holding period assuming that P/L are normally distributed, for specified confidence level and 
#'  holding period.
#' 
#' @param ... The input arguments contain either return data or else mean and                                                  
#'  standard deviation data. Accordingly, number of input arguments is either 3 
#'  or 4. In case there 3 input arguments, the mean and standard deviation of 
#'  data is computed from return data. See examples for details.
#' 
#' returns Vector of daily geometric return data
#' 
#'  mu Mean of daily geometric return data
#' 
#'  sigma Standard deviation of daily geometric return data
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
#'    data <- rnorm(5, .07, .03)
#'    NormalVaRPlot3D(returns = data, cl = seq(.9,.99,.01), hp = 1:100)
#'    
#'    # Computes VaR against confidence level given mean and standard deviation of return data
#'    NormalVaRPlot3D(mu = .012, sigma = .03, cl = seq(.9,.99,.01), hp = 1:100)
#'
#'
#' @export
NormalVaRPlot3D <- function(...){
  # Determine if there are three or four arguments, and ensure that arguments are read as intended
  if (nargs() < 3) {
    stop("Too few arguments")
  }
  if (nargs() > 4) {
    stop("Too many arguments")
  }
  args <- list(...)
  if (nargs() == 4) {
    mu <- args$mu
    cl <- args$cl
    sigma <- args$sigma
    hp <- args$hp
  }
  if (nargs() == 3) {
    mu <- mean(args$returns)
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
  
  # Check that cl is read as row vector
  if (cl.row > cl.col) {
    cl <- t(cl)
  }
  # Check that hp is read as column vector
  if (hp.col > hp.row) {
    hp <- t(hp)
  }
  
  # Check that inputs obey sign and value restrictions
  if (sigma < 0) {
    stop("Standard deviation must be non-negative")
  }
  if (max(cl) >= 1){
    stop("Confidence levels must be less than 1")
  }
  if (min(cl) <= 0){
    stop("Confidence levels must be greater than 0")
  }
  if (min(hp) <= 0){
    stop("Holding periods must be greater than 0")
  }
  
  # VaR estimation
  cl.row <- dim(cl)[1]
  cl.col <- dim(cl)[2]
  VaR <- - sigma[1,1] * sqrt(hp) %*% qnorm(1 - cl, 0, 1)  - mu[1,1] * hp %*% matrix(1,cl.row,cl.col) # VaR
  # Plotting
  persp(x=cl, y=hp, t(VaR), xlab = "Confidence Level", 
        ylab = "Holding Period", zlab = "VaR", border=NA,
        theta = -45, phi = 40, shade = .75, ltheta = 120, cex.axis=.85, cex.lab=.85,
        col = "lightgray", ticktype = "detailed", nticks = 5, 
        main = "Normal VaR against CL and HP")
  
}
