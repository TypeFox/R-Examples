#' @export
#' @title Rolling Median
#' @param x an \R numeric vector
#' @param n integer window size
#' @param increment integer shift to use when sliding the window to the next location
#' @description Fast, center-aligned rolling medians using C++/Rcpp.
#' Additional performance gains can be achieved by skipping \code{increment} values between calculations.
#' 
#' The \code{roll_median} function can be used to replace outliers detected by the \code{roll_hampel} function.  See example below.
#' @details The window size \code{n} is interpreted as the full window length.
#' Values within \code{n/2} of the beginning or end of \code{x} are set to \code{NA}.
#'   
#' Setting \code{increment} to a value greater than one will result in \code{NA}s for all skipped-over indices.
#' @return A vector of rolling median values of the same length as \code{x}.
#' @seealso \code{\link{roll_hampel}}
#' @examples
# Noisy sinusoid with outliers
#' a <- jitter(sin(0.1*seq(1e4)),amount=0.2)
#' indices <- sample(seq(1e4),20)
#' a[indices] <- a[indices]*10
#' 
#' # Outlier detection
#' b <- roll_hampel(a,10)
#' threshold <- 6
#' outliers <- which(b > threshold)
#' 
#' # Outlier replacement with median values
#' a_fixed <- a
#' a_fixed[outliers] <- roll_median(a,10)[outliers]
#' 
#' # Set up two plots
#' layout(matrix(seq(2)))
#' plot(a,type='l', col='gray60', main="Outliers detected")
#' points(outliers,a[outliers], col='red', lwd=2)
#' 
#' plot(a_fixed,type='l', col='gray60',
#'      main="Outliers replaced with rolling median")
#' points(outliers,a_fixed[outliers], col='red', lwd=2)

roll_median <- function( x, n=7, increment=1 ) {
  
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

    result <- .Call("seismicRoll_roll_median_numeric_vector", 
                    x, as.integer(n), as.integer(increment),
                    PACKAGE="seismicRoll")
    
    return (as.numeric(result))
    
  }

}

