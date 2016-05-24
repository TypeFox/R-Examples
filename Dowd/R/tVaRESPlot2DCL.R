#' Plots t VaR and ES against confidence level
#' 
#' Plots the VaR and ES of a portfolio against confidence level assuming that P/L
#' data are t distributed, for specified confidence level and 
#'  holding period.
#' 
#' @param ... The input arguments contain either return data or else mean and 
#'  standard deviation data. Accordingly, number of input arguments is either 4
#'  or 5. In case there are 4 input arguments, the mean and standard deviation of 
#'  data is computed from return data. See examples for details.
#' 
#'  returns Vector of daily geometric return data
#' 
#'  mu Mean of daily geometric return data
#' 
#'  sigma Standard deviation of daily geometric return data
#' 
#'  cl VaR confidence level and must be a vector
#' 
#'  hp VaR holding period and must be a scalar
#'  
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#'
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Plots VaR and ETL against confidene level given P/L data
#'    data <- runif(5, min = 0, max = .2)
#'    tVaRESPlot2DCL(returns = data, df = 7, cl = seq(.85,.99,.01), hp = 60)
#'    
#'    # Computes VaR against confidence level given mean and standard deviation of P/L data
#'    tVaRESPlot2DCL(mu = .012, sigma = .03, df = 7, cl = seq(.85,.99,.01), hp = 40)
#'
#'
#' @export
tVaRESPlot2DCL<- function(...){
  # Determine if there are four or five arguments, and ensure that arguments 
  # are read as intended
  if (nargs() < 4) {
    stop("Too few arguments")
  }
  if (nargs() > 5) {
    stop("Too many arguments")
  }
  args <- list(...)
  if (nargs() == 5) {
    mu <- args$mu
    cl <- args$cl
    sigma <- args$sigma
    hp <- args$hp
    df <- args$df
  }
  if (nargs() == 4) {
    mu <- mean(args$returns)
    cl <- args$cl
    sigma <- sd(args$returns)
    hp <- args$hp
    df <- args$df
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
    stop("Holding period must be greater than 0")
  }
  # VaR estimation  
  cl.row <- dim(cl)[1]
  cl.col <- dim(cl)[2]
  VaR <- (-sigma[1,1] * sqrt(t(hp)) %*% sqrt((df - 2) / df) %*% qt(1 - cl, df)) + (- mu[1,1] * t(hp) %*% matrix(1, cl.row, cl.col)) # VaR
  
  # ES estimation
  n <- 1000 # Number of slices into which tail is divided
  cl0 <- cl # Initial confidence level
  v <- VaR 
  delta.cl <- (1 - cl)/n # Increment to confidence level as each slice is taken
  for (i in 1:(n-1)) {
    cl <- cl0 + i * delta.cl # Revised cl
    v <- v + (-sigma[1,1] * sqrt(t(hp)) %*% sqrt((df - 2) / df) %*% qt(1 - cl, df)) + (- mu[1,1] * t(hp) %*% matrix(1, cl.row, cl.col))
  }
  v <- v/n # ES
    
  # Plotting
  ymin <- min(VaR, v)
  ymax <- max(VaR, v)
  xmin <- min(cl0)
  xmax <- max(cl0)
  plot(cl0, VaR, type = "l", xlim = c(xmin, xmax), ylim = c(ymin, ymax), xlab = "Confidence level", ylab = "VaR/ETL")
  par(new=TRUE)
  plot(cl0, v, type = "l", xlim = c(xmin, xmax), ylim = c(ymin, ymax), xlab = "Confidence level", ylab = "VaR/ETL")
  
  title("t VaR and ETL against confidence level")
  xmin <- min(cl0)+.3*(max(cl0)-min(cl0))
  text(xmin,max(VaR)-.1*(max(VaR)-min(VaR)),
       'Input parameters', cex=.75, font = 2)
  text(xmin,max(VaR)-.175*(max(VaR)-min(VaR)),
       paste('Daily mean L/P = ',round(mu[1,1],3)),cex=.75)
  text(xmin,max(VaR)-.25*(max(VaR)-min(VaR)),
       paste('Stdev. of daily L/P = ',round(sigma[1,1],3)),cex=.75)
  text(xmin,max(VaR)-.325*(max(VaR)-min(VaR)),
       paste('Degrees of freedom = ',df),cex=.75)
  text(xmin,max(VaR)-.4*(max(VaR)-min(VaR)),
       paste('Holding period = ',hp,'days'),cex=.75)
  # VaR and ETL labels
  text(max(cl0)-.4*(max(cl0)-min(cl0)),min(VaR)+.3*(max(VaR)-min(VaR)),'Upper line - ETL',cex=.75);
  text(max(cl0)-.4*(max(cl0)-min(cl0)),min(VaR)+.2*(max(VaR)-min(VaR)),'Lower line - VaR',cex=.75);

}
