#' @export
#' @title Outlier Detection with a Rolling Hampel Filter
#' @param x an \R numeric vector
#' @param n integer window size
#' @param thresholdMin minimum value for outlier detection
#' @param selectivity value between [0-1] used in determining outliers
#' @param increment integer shift to use when sliding the window to the next location
#' @description   A wrapper for the roll_hampel() function that applies a simple
#' algorithm for determining appropriate threshold settings based on 
#' statistics of the incoming data.
#' @details   The \code{thresholdMin} level is similar to a sigma value for normally distributed data.
#' Hampel filter values above 6.0 indicate a data value that is extremely unlikely
#' to be part of a normal distribution  (~ 1/500 million) and therefore very likely to be an outlier. By
#' choosing a relatively large value for \code{thresholdMin} we make it less likely that we
#' will generate false positives.
#' 
#' The \code{selectivity} is a value between 0 and 1 and is used to generate an appropriate 
#' threshold for outlier detection based on the statistics of the incoming data. A lower value
#' for \code{selectivity} will result in more outliers while a value closer to 1.0 will result in 
#' fewer.
#' 
#' The \code{thresholdMin} and \code{selectivity} parameters work like squelch and 
#' volume on a CB radio: \code{thresholdMin} sets a noise threshold below which you don't want anything
#' while \code{selectivity} increases the number of points defined as outliers. Of course \code{n},
#' the windowSize, is important as well.

#' @note For high resolution seismic data available from the 
#' \href{http://www.iris.edu/dms/nodes/dmc/}{IRIS DMC}, \code{B..} or \code{V..} channels, the default value
#' of \code{thresholdMin=6.0} seems to work well. For 1 Hz, \code{L..} channels a value
#' of \code{thresholdMin=12.0} seems more appropriate.
#' The user is advised to do some testing before applying this filter to low resolution
#' seismic data.
#' 
#' \bold{The default value of \code{increment=1} should not be changed.} Outliers are defined
#' as individual points that stand apart from their neighbors. Applying the Hampel filter to
#' every other point by using \code{increment} > 1 will invariably miss some of the outliers.
#' @return A vector of indices associated with outliers in the incoming data \code{x}.
#' @seealso \code{\link{roll_hampel}}
#' @examples
#' # Noisy sinusoid with outliers
#' a <- jitter(sin(0.1*seq(1e4)),amount=0.2)
#' indices <- sample(seq(1e4),20)
#' a[indices] <- a[indices]*10
#' 
#' # Outlier detection should identify many of these altered indices
#' sort(indices)
#' findOutliers(a)

findOutliers <- function(x, n=41, thresholdMin=6, selectivity=0.4, increment=1) {
  
  h <- roll_hampel(x, n, increment)

  maxH <- max(h, na.rm=TRUE)
  
  if (maxH < thresholdMin) {
    # If no values cross thresholdMin, return an empty vector
    
    return(which(h > maxH)) # integer vector of length zero
    
  } else {
    # Some values cross thresholdMin.
    # Set up a new threshold based on maxH and selectivity
    
    outliers <- which(h > maxH * selectivity)
    return(outliers)
  }
  
}

