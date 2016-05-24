#' @title Normalize a Numeric Vector
#' @description Normalize a numeric vector to have a mean of zero and a standard deviation of one.\
#' 
#' @param x a numeric vector.
#' @return a numeric vector of the same length as \code{x}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' x <- runif(20, 50, 110)
#' x.norm <- normalize(x)
#' mean(x)
#' mean(x.norm)
#' sd(x)
#' sd(x.norm)
#' 
#' @importFrom stats sd
#' @export
#' 
normalize <- function(x) {
  if(!is.numeric(x) & !is.vector(x)) stop("'x' must be a numeric vector")
  x.mean <- mean(x, na.rm = TRUE)
  x.sd <- sd(x, na.rm = TRUE)
  (x - x.mean) / x.sd
}
