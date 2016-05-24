#' @title Percentile of historical simulation VaR distribution function
#'
#' @description Estimates percentiles of historical simulation VaR distribution 
#' function, using theory of order statistics, for specified confidence level.
#'
#' @param Ra Vector of daily P/L data
#' @param perc Desired percentile and is scalar
#' @param cl VaR confidence level and is scalar
#' @return Value of percentile of VaR distribution function
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Estimates Percentiles for random standard normal returns and given perc 
#'    # and cl
#'    Ra <- rnorm(100)
#'    HSVaRDFPerc(Ra, .75, .95)
#'
#' @export
HSVaRDFPerc <- function(Ra, perc, cl){
  
  # Determine if there are three arguments, and ensure that arguments are read as intended
  if (nargs() < 3) {
    stop("Too few arguments.")
  }
  if (nargs() > 3) {
    stop("Too many arguments")
  }
  if (nargs() == 3) {
    profit.loss <- as.vector(Ra)
    data <- sort(profit.loss)
    n <- length(data)
  }
  
  # Check that inputs obey sign and value restrictions
  if (n < 0) {
    stop("Number of observations must be greater than zero.")
  }
  if (perc <= 0) {
    stop("Chosen percentile must be positive.")
  }
  if (perc > 1) {
    stop("Chosen percentile must not exceed 1")
  }
  if (cl >= 1) {
    stop("Confidence level must be less than 1.")
  }
  if (cl <= 0) {
    stop("Confidence level must positive.")
  }
  
  # Derive order statistics and ensure it is an integer
  w <- n * cl # Derive rth order statistics
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
  i <- round(n * x)
  y <- data[i] # Value of percentile of VaR distribution function
  return(y)
  
} 