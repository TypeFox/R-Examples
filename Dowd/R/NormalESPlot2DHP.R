#' Plots normal ES against holding period
#' 
#' Plots the ES of a portfolio against holding period assuming that P/L distribution is 
#' normally distributed, for specified confidence level and holding period.
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
#'  cl ES confidence level and must be a scalar
#'  
#'  hp ES holding period and must be a vector
#'  
#'  
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#'
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Computes ES given geometric return data
#'    data <- runif(5, min = 0, max = .2)
#'    NormalESPlot2DHP(returns = data, cl = .95, hp = 60:90)
#'    
#'    # Computes v given mean and standard deviation of return data
#'    NormalESPlot2DHP(mu = .012, sigma = .03, cl = .99, hp = 40:80)
#'
#' @export
NormalESPlot2DHP <- function(...){
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
    stop("Confidence levels must be less than 1")
  }
  if (min(cl) <= 0){
    stop("Confidence levels must be greater than 0")
  }
  if (min(hp) <= 0){
    stop("Holding Period must be greater than 0")
  }
  # VaR estimation
  
  cl.row <- dim(cl)[1]
  cl.cl <- dim(cl)[2]
  VaR <- - sigma[1,1] * sqrt(t(hp)) * qnorm(1 - cl[1,1], 0, 1)
                          + mu[1,1] * t(hp) %*% matrix(1, cl.row, cl.col) # VaR
  
  # ES etimation
  es <- NormalES(mu = mu[1,1], sigma = sigma[1,1], cl = cl, hp = hp)
  
  # Plotting
  plot(hp, es, type = "l", xlab = "Holding Period", ylab = "ES")
  title("Normal ES against holding period")
  cl.label <- 100*cl
  xmin <-min(hp)+.25*(max(hp)-min(hp))
  text(xmin,max(es)-.1*(max(es)-min(es)),
       'Input parameters', cex=.75, font = 2)
  text(xmin,max(es)-.15*(max(es)-min(es)),
       paste('Daily mean L/P = ',mu[1,1]),cex=.75)
  text(xmin,max(es)-.2*(max(es)-min(es)),
       paste('Stdev. of daily L/P = ',sigma[1,1]),cex=.75)
  text(xmin,max(es)-.25*(max(es)-min(es)),
       paste('Confidence level = ',cl.label,'%'),cex=.75)
}
