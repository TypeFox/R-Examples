#' Value at Risk of a portfolio using Historical Estimator
#'
#' Estimates the Value at Risk (VaR) using historical estimator 
#' approach for the specified range of confidence levels and the holding
#' period implies by data frequency.
#'
#' @param Ra Vector corresponding to profit and loss distribution
#' @param Rb Scalar corresponding to VaR confidence levels.
#' @return Value at Risk of the portfolio
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' Jorion, P. Value at Risk: The New Benchmark for Managing Financial Risk. 
#' McGraw-Hill, 2006
#' 
#' Cont, R., Deguest, R. and Scandolo, G. Robustness and sensitivity analysis
#' of risk measurement procedures. Quantitative Finance, 10(6), 2010, 593-606.
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
#'    # To be added
#'    a <- rnorm(1000) # Payoffs of random portfolio
#'    HSVaR(a, .95)
#'
#' @export
HSVaR <- function(Ra, Rb){
  
  # Determine if there are two arguments and ensure that they are read as
  # intended
  if (nargs() < 2) {
    stop("Too few arguments")
  }
  
  if (nargs() > 2) {
    stop("Too many arguments")
  }
  
  if (nargs() == 2) {
    profit.loss.data <- as.vector(Ra)
    cl <- as.vector(Rb)
    unsorted.loss.data <- -profit.loss.data # Derives L/P data from input P/L
    losses.data <- sort(unsorted.loss.data) # Puts losses in ascending order
    n <- length(losses.data)
  }
  
  # Check that inputs have correct dimensions
  if (!is.vector(cl)) {
    stop("Confidence level must be a vector")
  }
  
  # Check that inputs obey sign and value restrictions
  if (max(cl) >= 1) {
    stop("Confidence level must be less than 1.")
  }
  if (max(cl) <= 0) {
    stop("Confidence level must be positive")
  }
  
  i <- 1:length(cl)
  index <- cl*n # This putative index value may or may not be an integer
  
  # If index value is an integer, VaR follows immediately
  y <- double(length(i))
  if (index-round(index) == 0){
    y[i] <- losses.data[index]
  }
  
  # If index not an integer, take VaR as linear interpolation of loss
  # observations just above and below "true" VaR
  
  if (index-round(index)!=0){
    # Deal with loss observation just above VaR
    upper.index <- ceiling(index)
    upper.var <- losses.data[upper.index] # Loss observation just above VaR or upper VaR
    
    # Deal with loss observation just below VaR
    lower.index <- floor(index)
    lower.var <- losses.data[lower.index] # Loss observation just below VaR or lower VaR
    
    # If lower and upper indices are the same ,VaR is upper VaR
    if (upper.index==lower.index){
      y <- upper.var
    }
    
    # If lower and upper indices different, VaR is weighted average of upper
    # and lower VaRs
    if (upper.index!=lower.index){
      # Weights attached to upper and lower VaRs
      lower.weight <- (upper.index-index)/(upper.index-lower.index) # weight on upper.var
      upper.weight <- (index-lower.index)/(upper.index-lower.index) # weight on upper_var
      # Finally, the weighted, VaR as a linear interpolation of upper and lower VaRs
      
      y <- lower.weight * lower.var + upper.weight * upper.var
      
    }
    
  }
    return(y)
}
