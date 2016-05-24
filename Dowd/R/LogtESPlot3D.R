#' Plots log-t ES against confidence level and holding period
#' 
#' Plots the ES of a portfolio against confidence level and holding period assuming that geometric
#'  returns are Student-t distributed, for specified confidence level and 
#'  holding period.
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
#'  cl VaR confidence level and must be a vector
#' 
#'  hp VaR holding period and must be a vector
#'  
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#'
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Plots ES against confidene level given geometric return data
#'    data <- rnorm(5, .09, .03)
#'    LogtESPlot3D(returns = data, investment = 5, df = 6, cl = seq(.9,.99,.01), hp = 1:100)
#'    
#'    # Computes ES against confidence level given mean and standard deviation of return data
#'    LogtESPlot3D(mu = .012, sigma = .03, investment = 5, df = 6, cl = seq(.9,.99,.01), hp = 1:100)
#'
#'
#' @export
LogtESPlot3D <- function(...){
  # Determine if there are four or five arguments, and ensure that arguments are read as intended
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
    stop("Confidence level(s) must be less than 1")
  }
  if (min(cl) <= 0){
    stop("Confidence level(s) must be greater than 0")
  }
  if (min(hp) <= 0){
    stop("Confidence level(s) must be greater than 0")
  }
  
  # VaR estimation
  cl.row <- dim(cl)[1]
  cl.col <- dim(cl)[2]
  hp.row <- dim(hp)[1]
  hp.col <- dim(hp)[2]
  
  VaR <- investment - exp( ((df-2)/df) * sigma[1,1] * sqrt(hp) %*% qt(1 - cl, df)  + mu[1,1] * hp %*% matrix(1,cl.row,cl.col) + log(investment)) # VaR
  
  # ES etimation
  n <- 1000 # Number of slices into which tail is divided
  cl0 <- cl # Initial confidence level
  delta.cl <- (1 - cl) / n # Increment to confidence level as each slice is taken
  v <- VaR
  for (i in 1:(n-1)) {
    cl <- cl0 + i * delta.cl # Revised cl
    v <- v + investment - exp( ((df-2)/df) * sigma[1,1] * sqrt(hp) %*% qt(1 - cl, df)  + mu[1,1] * hp %*% matrix(1,cl.row,cl.col) + log(investment))
  }
  v <- v/n
  
  # Plotting
  persp(x=cl, y=hp, t(v), xlab = "Confidence Level", 
        ylab = "Holding Period", zlab = "ES", border=NA,
        theta = -45, phi = 35, shade = .75, ltheta = 90, cex.axis=.85, cex.lab=.85,
        col = "lightgray", ticktype = "detailed", nticks = 5,
        main = "Log-t ES against CL and HP")
  
}
