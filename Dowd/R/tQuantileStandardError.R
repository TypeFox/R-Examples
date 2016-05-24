#' Standard error of t quantile estimate
#'
#' Estimates standard error of t quantile estimate
#'  
#' @param prob Tail probability. Can be a vector or scalar
#' @param n Sample size
#' @param mu Mean of the normal distribution
#' @param sigma Standard deviation of the distribution
#' @param df Number of degrees of freedom
#' @param bin.size Bin size. It is optional parameter with default value 1
#' @return Vector or scalar 
#' depending on whether the probability is a vector
#' or scalar
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#'
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Estimates standard error of normal quantile estimate
#'    tQuantileStandardError(.8, 100, 0, .5, 5, 3)
#'
#' @export
tQuantileStandardError <- function(prob, n, mu, sigma, df, bin.size){
  # Check that inputs obey sign and value restrictions
  if (prob < 0|prob>1) {
    stop("Probability must be nonnegative and no greater than 1")
  }
  if (n <= 0){
    stop("Sample size must be positive")
  }
  if (bin.size <= 0){
    stop("Bin size must be greater than 0")
  }
  # Determination of frequency 
  x <- mu + sigma * qt(prob, df)
  z <- (x - mu)/sigma
  freq <- pt((x + .5 * bin.size - mu) / sigma, df) - pt((x - 0.5 * bin.size - mu) / sigma, df)
  #  Standard error estimation
  y <- prob * (1 - prob) / (n * freq ^ 2) # Standard Error
  return(y)
}