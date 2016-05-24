#' @export
#' @title Rolling Hampel Filter for Outlier Detection
#' @param x an \R numeric vector
#' @param n integer window size
#' @param increment integer shift to use when sliding the window to the next location
#' @description Fast, center-aligned hampel filter using C++/Rcpp.
#' The Hampel filter is a robust outlier detector using Median Absolute Deviation (MAD).
#' Additional performance gains can be achieved by skipping \code{increment} values between calculations. 
#' @details   \emph{Unlike} the version in the \pkg{pracma} package, this version does not return the
#' corrected timeseries. Instead, it returns a vector of values that can be tested against
#' different threshold values. Higher values in the return are associated with a higher likelihood
#' that the associated point is an outlier when compared with its neighbors. Outliers can
#' be picked out by comparing the return values against some threshold as seen in the example.
#' 
#' Also \emph{unlike} the \pkg{pracma} version, \code{n} is interpreted as the full window length
#' and will be increased by one if necessary to have a window representing an odd number of indices.
#' @note A pure \R version of the filter is found in the \pkg{pracma} package.
#' @return A vector of values of the same length as \code{x}.
#' @seealso \code{\link{roll_median}}
#' @examples
#' a <- sin(0.1*seq(100))
#' a[20] <- 50
#'
#' b <- roll_hampel(a,10)
#' threshold <- 6
#' which(b > threshold)
#'
#' \dontrun{
#'   require(microbenchmark)
#'   require(pracma)
#'   
#'   microbenchmark(hampel(a,10), roll_hampel(a,10), times=10)
#'   
#'   #  Unit: microseconds
#'   #                 expr      min       lq    median       uq       max neval
#'   #        hampel(a, 10) 7610.688 7784.374 8037.4035 9453.928 16176.535    10
#'   #   roll_hampel(a, 10)   36.530   37.443   58.7165   65.418    90.403    10
#' }

roll_hampel <- function( x, n, increment=1 ) {
  
  if ( !is.vector(x) ) {

    stop("the x supplied is not a vector")
    
  } else {
    
    if ( n > length(x) ) {
      stop("n cannot be greater than length(x).")
    }

    # Avoid infinite loop
    if ( increment < 1 ) {
      stop("increment must be >= 1.")
    }

    result <- .Call( "seismicRoll_roll_hampel_numeric_vector", 
                     x, as.integer(n), as.integer(increment), 
                     PACKAGE="seismicRoll")
    
    return (as.numeric(result))
    
  }
  
}

