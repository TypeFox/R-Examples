#' @title Unbiased Estimate of Diversity
#' @description Calculate unbiased estimate of diversity for a vector of items
#' 
#' @param x charcter or numeric vector or factor 
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @importFrom stats na.omit
#' @export
#' 
diversity <- function(x) {
  if(!(is.vector(x) | is.factor(x))) {
    stop("'x' must be a character or numeric vector, or a factor")
  }
  x <- na.omit(x)
  x.freq <- prop.table(table(x))
  n <- length(x)    
  n * (1 - sum(x.freq ^ 2)) / (n - 1)
}