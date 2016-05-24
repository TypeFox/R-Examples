#' Plots normal ES against confidence level
#' 
#' Plots the ES of a portfolio against confidence level assuming that P/L are 
#' normally distributed, for specified confidence level and holding period.
#' 
#' @param ... The input arguments contain either return data or else mean and 
#'  standard deviation data. Accordingly, number of input arguments is either 3 
#'  or 4. In case there 3 input arguments, the mean and standard deviation of 
#'  data is computed from return data. See examples for details.
#'  
#'  returns Vector of daily geometric return data
#'  
#'  mu Mean of daily geometric return data
#'  
#'  sigma Standard deviation of daily geometric return data
#'  
#'  cl ES confidence level and must be a vector
#'  
#'  hp ES holding period and must be a scalar
#'  
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#'
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Plots ES against confidence level
#'    data <- runif(5, min = 0, max = .2)
#'    NormalESPlot2DCL(returns = data, cl = seq(.9,.99,.01), hp = 60)
#'    
#'    # Plots ES against confidence level
#'    NormalESPlot2DCL(mu = .012, sigma = .03, cl = seq(.9,.99,.01), hp = 40)
#'
#' @export
NormalESPlot2DCL <- function(...){
  # Determine if there are four or five arguments, and ensure that arguments are read as intended
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
  if (max(hp.row, hp.col) > 1) {
    stop("Holding period must be a scalar")
  }
  
  # Check that hp is read as row vector
  if (cl.row > cl.col) {
    cl <- t(cl)
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
    stop("Holding period must be greater than 0")
  }

  # ES estimation
  v <- NormalES(mu = mu[1,1], sigma = sigma[1,1], cl = as.vector(cl), hp = hp[1,1])
  
  # Plotting
  plot(cl, v, type = "l", xlab = "Holding Period", ylab = "ES")
  title("Normal ES against holding period")
  xmin <-min(cl)+.25*(max(cl)-min(cl))
  text(xmin,max(v)-.1*(max(v)-min(v)),
       'Input parameters', cex=.75, font = 2)
  text(xmin,max(v)-.15*(max(v)-min(v)),
       paste('Daily mean geometric return = ',round(mu[1,1],3)),cex=.75)
  text(xmin,max(v)-.2*(max(v)-min(v)),
       paste('Stdev. of daily geometric returns = ',round(sigma[1,1],3)),cex=.75)
  text(xmin,max(v)-.25*(max(v)-min(v)),
       paste('Holding Period = ',hp),cex=.75)
}
