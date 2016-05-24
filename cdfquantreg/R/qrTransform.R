#' @title Transform Values into (0, 1) Interval
#' @aliases scaleTR
#' @description \code{scaleTR} is function that rescales values of a variable into the (0, 1) interval.
#' @param y A numeric vector, or a variable in a dataframe. 
#' @param high The highest possible value of that variable. The value should be equal or greater than the maximum value of y. If not supplied, the maximum value of y will be used. 
#' @param low The lowest possible value of that variable. The value should be equal or smaller than the minimum value of y. If not supplied, the minimum value of y will be used.  
#' @param data A dataframe that contains the variable y.  
#' @param N A integer, normally is the sample size or the number of values. If not supplied, the length of y will be used. 
#' @param scale A compressing parameter that determines the extend to which the boundary values are going to be pushed away from the boundray. See details.
#' 
#' @details \code{scaleTR} used the method suggested by Smithson and Verkuilen (2006) and applies linear transformation to values into the open interval (0, 1). It first transform the values from their original scale by taking \eqn{y' = (y - a)/(b-a)}, where \code{a} is the lowest possible value of that variable and  \code{b} is the highest possible value of that variable. Next, it compresses the range to avoid zeros and ones by taking \eqn{y" = (y'(N - 1) + c)/N}, where \code{N} is the sample size and \code{c} is the compressing parameter. The smaller value \code{c} is, the boundray values would be more approaching zeros and ones, and have greater impact on the estimation of the dispersion parameters in the cdf quantile model. 
#' 
#' @seealso \code{\link{cdfquantreg}}
#' @export 
#' 
#' @examples
#' y <- rnorm(20, 0, 1)
#' ynew <- scaleTR(y)
#' 

scaleTR <- function(y, high = NULL, low = NULL, data = NULL,  
                       N = NULL, scale = 0.5) {
  
  if (!is.null(data)) {
    yname <- deparse(substitute(y))
    y <- data[, yname]
  }

  if (is.null(high)) {
    high <- max(y)
  }
  
  if (is.null(low)) {
    low <- min(y)
  }
  
  if (is.null(N)) {
    N <- length(y)
  }
  
  y1 <- (y - low)/(high-low)
  
  y2 <- (y1*(N - 1) + scale)/N
  
  if (is.null(data)) {
    return(y2)
  } else {
   names <- paste(yname,'old', sep ='')
   data[, names] <- y
   data[, yname] <- y2
   return(data)
  }

}

