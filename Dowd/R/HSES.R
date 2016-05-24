#' Expected Shortfall of a portfolio using Historical Estimator
#'
#' Estimates the Expected Shortfall (aka. Average Value at Risk or Conditional 
#' Value at Risk) using historical estimator approach for the specified
#' confidence level and the holding period implies by data frequency.
#'
#' @param Ra Vector corresponding to profit and loss distribution
#' @param cl Number between 0 and 1 corresponding to confidence level
#' @return Expected Shortfall of the portfolio
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' Cont, R., Deguest, R. and Scandolo, G. Robustness and sensitivity analysis
#' of risk measurement procedures. Quantitative Finance, 10(6), 2010, 593-606.
#' 
#' Acerbi, C. and Tasche, D. On the coherence of Expected Shortfall. Journal
#' of Banking and Finance, 26(7), 2002, 1487-1503
#' 
#' Artzner, P., Delbaen, F., Eber, J.M. and Heath, D. Coherent Risk Measures 
#' of Risk. Mathematical Finance 9(3), 1999, 203.
#' 
#' Foellmer, H. and Scheid, A. Stochastic Finance: An Introduction in Discrete 
#' Time. De Gryuter, 2011.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Computes Historical Expected Shortfall for a given profit/loss
#'    # distribution and confidence level
#'    a <- rnorm(100) # generate a random profit/loss vector
#'    HSES(a, 0.95)
#'
#' @export
HSES <- function(Ra, cl){
  
  if (nargs() < 2) {
    stop("Too few arguments")
  }
  
  if (nargs() > 2) {
    stop("Too many arguments")
  }
  
  if (nargs() == 2) {
    profit.loss.data <- as.vector(Ra)
    unsorted.loss.data <- -profit.loss.data # Derives L/P data from input P/L
    losses.data <- sort(unsorted.loss.data) # Puts losses in ascending order
    n <- length(losses.data)
  }
  
  # Check that inputs have correct dimensions
  if (length(cl) != 1) {
    stop('Confidence level must be scalar (length-1 vector in R)')
  }
  
  # Check that inputs obey sign and value restrictions
  if (cl >= 1) {
    stop("Confidence level must be less than 1.")
  }
  if (cl <= 0) {
    stop("Confidence level must be positive")
  }
  
  # VaR and ES estimation
  index <- n*cl # This putative index value may or may not be an integer
  
  # Each case needs to be considered in turn
  # If index value is an integegr, VaR follows immediately and then we
  # estimate ES
  if (index-round(index)==0){
    VaR <- losses.data[index] # Historical Value at Risk
    k <- which(VaR <= losses.data) # Finds indices of tail loss data
    tail.losses <- losses.data[k] # Creates data set of tail loss observations
    es <- mean(tail.losses) # Expected Shortfall
    y <- es    
  }
  
  # If index not an integer, take VaR as linear interpolation of loss
  # observationsjust above and "below" true VaR and take Expected Shortfall
  # as linear interpolation of corresponding upper and lower Expected Shortfall
  if (index-round(index) != 0){
    # Deal with loss
    upper.index <- ceiling(index)
    upper.VaR <- losses.data[upper.index] # Upper VaR
    upper.k <- which(upper.VaR<=losses.data) # Finds indices of upper tail loss data
    upper.tail.losses <- losses.data[upper.k] # Creates data set of upper tail loss obs.
	upper.es <- mean(upper.tail.losses) # Upper ES
	# Deal with loss observation just below VaR to derive lower ES
	lower.index <- ceiling(index)
	lower.VaR <- losses.data[lower.index] # Lower VaR
	lower.k <- which(lower.VaR <= losses.data) # Finds indices of lower tail loss data 
	lower.tail.losses <- losses.data[lower.k] # Creates data set of lower tail loss obs.
	lower.es <- mean(lower.tail.losses)# Lower ES
	
    lower.es <- mean(lower.tail.losses) # Lower Expected Shortfall (ES)
    # If lower and upper indices are the same, ES is upper ES
    if (upper.index == lower.index){
      y <- upper.es
    }
    # If lower and upper indices are different, ES is weighted average of
    # upper and lower ESs
    if (upper.index!=lower.index) {
      # Weights attached to upper and lower ESs
      lower.weight <- (upper.index-index)/(upper.index-lower.index) # weight on upper_var
      upper.weight <- (index-lower.index)/(upper.index-lower.index) # weight on upper_var
      # Finally, the weighted, ES as a linear interpolation of upper and lower
      # ESs
      y <- lower.weight*lower.es+upper.weight*upper.es
    }
  }
  return(y)
}
