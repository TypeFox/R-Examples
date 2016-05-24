#' Plots log normal VaR against holding period
#' 
#' Plots the VaR of a portfolio against holding period assuming that geometric returns are 
#' normal distributed, for specified confidence level and holding period.
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
#'  investment Size of investment
#' 
#'  cl VaR confidence level and must be a scalar
#' 
#'  hp VaR holding period and must be a vector
#'  
#'  @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#'
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Computes VaR given geometric return data
#'    data <- runif(5, min = 0, max = .2)
#'    LogNormalVaRPlot2DHP(returns = data, investment = 5, cl = .95, hp = 60:90)
#'    
#'    # Computes VaR given mean and standard deviation of return data
#'    LogNormalVaRPlot2DHP(mu = .012, sigma = .03, investment = 5, cl = .99, hp = 40:80)
#'
#'
#' @export
LogNormalVaRPlot2DHP <- function(...){
  # Determine if there are four or five arguments, and ensure that arguments are read as intended
  if (nargs() < 4) {
    stop("Too few arguments")
  }
  if (nargs() > 5) {
    stop("Too many arguments")
  }
  args <- list(...)
  if (nargs() == 5) {
    mu <- args$mu
    investment <- args$investment
    cl <- args$cl
    sigma <- args$sigma
    hp <- args$hp
  }
  if (nargs() == 4) {
    mu <- mean(args$returns)
    investment <- args$investment
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
  if (max(cl.row, cl.col) > 1) {
    stop("Confidence level must be a scalar")
  }
  hp <- as.matrix(hp)
  hp.row <- dim(hp)[1]
  hp.col <- dim(hp)[2]
  if (min(hp.row, hp.col) > 1) {
    stop("Holding period must be a vector")
  }
  
  # Check that hp is read as row vector
  if (hp.row > hp.col) {
    hp <- t(hp)
  }
  
  # Check that inputs obey sign and value restrictions
  if (sigma < 0) {
    stop("Standard deviation must be non-negative")
  }
  if (max(cl) >= 1){
    stop("Confidence level must be less than 1")
  }
  if (min(cl) <= 0){
    stop("Confidence level must be greater than 0")
  }
  if (min(hp) <= 0){
    stop("Holding periods must be greater than 0")
  }
  # VaR estimation  
  cl.row <- dim(cl)[1]
  cl.col <- dim(cl)[2]
  VaR <- investment - exp(sigma[1,1] * sqrt(t(hp)) * qnorm(1 - cl[1,1], 0, 1)
                        + mu[1,1] * t(hp) %*% matrix(1, cl.row, cl.col) + log(investment)) # VaR
  # Plotting
  plot(hp, VaR, type = "l", xlab = "Holding Period", ylab = "VaR")
  cl.label <- cl * 100
  title("Log Normal VaR against holding period")
  xmin <-min(hp)+.25*(max(hp)-min(hp))
  text(xmin,max(VaR)-.1*(max(VaR)-min(VaR)),
       'Input parameters', cex=.75, font = 2)
  text(xmin,max(VaR)-.15*(max(VaR)-min(VaR)),
       paste('Daily mean geometric return = ',mu[1,1]),cex=.75)
  text(xmin,max(VaR)-.2*(max(VaR)-min(VaR)),
       paste('Stdev. of daily geometric returns = ',sigma[1,1]),cex=.75)
  text(xmin,max(VaR)-.25*(max(VaR)-min(VaR)),
       paste('Investment size = ',investment),cex=.75)
  text(xmin,max(VaR)-.3*(max(VaR)-min(VaR)),
       paste('Confidence level = ',cl.label,'%'),cex=.75)
}
