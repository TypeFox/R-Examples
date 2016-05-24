#' # private functions
#' 

.optimalNumberOfBinsForHist <- function(x) {
  # Compute the optimal number of bins for a histogram.
  
  x.range.length <- diff(range(x))
  sd.x <- stats::sd(x)
  num.samples <- length(x)
  
  num.bins <- ceiling(x.range.length/(3.96 * sd.x * num.samples^(-1/3)))
  # num.bins <- nclass.FD(x)
  return(num.bins)
} 

