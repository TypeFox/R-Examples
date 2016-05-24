#' Figure of t- VaR and pdf against L/P
#'
#' Gives figure showing the VaR and probability distribution function against L/P
#'  of a portfolio assuming P/L are normally distributed, for specified 
#'  confidence level and holding period.
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
#'  df Number of degrees of freedom
#' 
#'  cl VaR confidence level and should be scalar
#' 
#'  hp VaR holding period in days and should be scalar
#'  
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#'
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Plots normal VaR and pdf against L/P data for given returns data
#'    data <- runif(5, min = 0, max = .2)
#'    tVaRFigure(returns = data, df = 7, cl = .95, hp = 90)
#'    
#'    # Plots normal VaR and pdf against L/P data with given parameters
#'    tVaRFigure(mu = .012, sigma = .03, df=7, cl = .95, hp = 90)
#'
#'
#' @export
tVaRFigure <- function(...){
  # Determine if there are four or five arguments and ensure that arguments are
  # read as intended
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
    stop("Confidence level must be a scalar or a vector")
  }
  hp <- as.matrix(hp)
  hp.row <- dim(hp)[1]
  hp.col <- dim(hp)[2]
  if (min(hp.row, hp.col) > 1) {
    stop("Holding period must be a scalar or a vector")
  }
  df <- as.matrix(df)
  df.row <- dim(df)[1]
  df.col <- dim(df)[2]
  if (max(df.row, df.col) > 1) {
    stop("Number of degrees of freedom must be a scalar")
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
  if (df < 3) {
    stop("Number of degrees of freedom must be at least 3 for first two moments
         of distribution to be defined")
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
  
  # Message to indicate how matrix of results is to be interpreted, if cl and hp
  # both vary and results are given in matrix form
  if (max(cl.row, cl.col) > 1 & max(hp.row, hp.col) > 1) {
    print('VaR results with confidence level varying across row and holding 
          period down column')
  }
  
  # VaR estimation
  cl.row <- dim(cl)[1]
  cl.col <- dim(cl)[2]
  VaR <- - sigma[1,1] * sqrt(hp) %*% qnorm(1 - cl, 0, 1)  - mu[1,1] * hp %*% matrix(1,cl.row,cl.col) # VaR
  
  # Plotting
  x.min <- -mu - 5 * sigma
  x.max <- -mu + 5 * sigma
  delta <- (x.max-x.min) / 100
  x <- seq(x.min, x.max, delta)
  p <- dt((x-mu)/sigma, df)
  plot(x, p, type = "l", xlim = c(x.min, x.max), ylim = c(0, max(p)*1.1), 
       xlab = "Loss (+) / Profit (-)", ylab = "Probability", main = "t- VaR")
  u <- c(VaR, VaR)
  v <- c(0, .6*max(p))
  lines(0,0,2,.6,type="l")
  lines(u, v, type = "l", col = "blue")
  cl.for.label <- 100*cl
  text(1,.95*max(p), pos = 1, 'Input parameters', cex=.75, font = 2)
  text(1, .875*max(p),pos = 1, paste('Daily mean L/P = ', round(mu,2)), cex=.75)
  text(1, .8*max(p),pos = 1, paste('St. dev. of daily L/P = ',round(sigma,2)), cex=.75)
  text(1, .725*max(p),pos = 1, paste('Holding period = ', hp,' day(s)'), cex=.75)
  text(VaR, .7*max(p),pos = 2, paste('VaR at ', cl.for.label,'% CL'), cex=.75)
  text(VaR, .64 * max(p),pos = 2, paste('= ',VaR), cex=.75)
}
