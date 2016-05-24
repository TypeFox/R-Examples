#' Harmonic mean
#'
#' Calculate the harmonic mean of a numeric vector
#' (will return NA if there are any negative numbers in the vector)
#' @param x numeric vector
#' @param na.rm logical remove NAs prior or calculation
#' @export
#' @export
#' @return harmonic mean of vector
#' @examples
#' 
#' data(nancycats)
#' pop.sizes <- table(pop(nancycats))
#' harmonic_mean(pop.sizes)


harmonic_mean <- function(x, na.rm=TRUE){
  if(na.rm){
      x <- x[!is.na(x)]
  } else {
      if(any(is.na(x))) return(NA)
  }
  if(any(x < 0)){
     return(NA)
  }
  1/mean(1/x)
}
