#' @title Harmonic Mean
#' @description Calculate the harmonic mean of a set of numbers.
#' 
#' @param x a numeric vector.
#' @param na.rm a logical value indicating whether NA values should be stripped 
#'   before the computation proceeds.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'  
#' @examples
#' x <- rlnorm(100)
#' mean(x)
#' median(x)
#' harmonic.mean(x)
#' 
#' @importFrom stats na.omit var
#' @export
#' 
harmonic.mean <- function(x, na.rm = FALSE) {
  #
  # Calculate harmonic mean of vector x
  # If any x < 0, then approximation used
  #
  #   9/15/2010
  
  if(na.rm) x <- na.omit(x)
  if(length(x) == 0) return(NA)
  hm <- if(all(x > 0)) {
    length(x) / sum(1 / x)
  } else {
    warning("some values are <= 0, using approximation")
    inv.mean.x <- 1 / mean(x)
    var.x <- var(x)
    1 / (inv.mean.x + var.x * inv.mean.x ^ 3)
  }
  ifelse(is.nan(hm), NA, hm)
}
