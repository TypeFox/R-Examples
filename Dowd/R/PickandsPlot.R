#' Pickand Estimator - Tail Sample Size Plot
#'
#' Displays a plot of the Pickands Estimator against Tail Sample Size.
#'
#' @param Ra The data set
#' @param maximum.tail.size maximum tail size and should be greater than a
#' quarter of the sample size.
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Pickand - Sample Tail Size Plot for random standard normal data
#'    Ra <- rnorm(1000)
#'    PickandsPlot(Ra, 40)
#'    
#'
#' @export
PickandsPlot <- function(Ra, maximum.tail.size){
  
  data <- as.vector(Ra)
  x <- sort(data)
  n <- length(x)
  # Derivation of Pickands estimators and tail size series
  pe <- double(maximum.tail.size)
  for(k in 2:maximum.tail.size){
    pe[k] <- PickandsEstimator(x, k)
  }
  # Plot of Pickands Estimator against tail size
  k <- 1:maximum.tail.size
  plot(k, pe, type="l", main = "Pickands Estimator against Tail Size",
       xlab = "Number of observations in tail (k)", 
       ylab = "Pickands Estimator")
  
}