#' @title Geometric Mean
#' @description Calculates the geometric mean of a vector.
#' 
#' @param x a numeric vector.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' x <- rlnorm(100)
#' mean(x)
#' median(x)
#' geometric.mean(x)
#' 
#' @importFrom stats na.omit
#' @export
#' 
geometric.mean <- function(x) {
  x <- na.omit(x)
  if(length(x) == 0) return(NA)
  prod(x) ^ (1 / length(x))
}