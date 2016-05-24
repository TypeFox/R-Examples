#' Standardize or Unstandarize the Column Range
#' 
#' Simple functions for rescaling a data matrix to a coded design and back. \code{stdrange} converts
#' the design in actual measurements into a coded design, while \code{ustdrange} reverses the process
#' (if the correct arguments are given).
#' @aliases ustdrange
#' @param x matrix containing the design, or an object coercible to a matrix.
#' @param mins vector of original values, one for each column, which should be recoded to the value -1;
#' or which have alreadty been recoded to -1. This and the next argument are both recycled if not of the correct length.
#' @param maxs vector of original values which should be recoded as 1, or which have already been recoded to 1.
#' @export
#' @author Pieter C. Schoonees

stdrange <- function(x, mins = apply(x, 2, min), maxs = apply(x, 2, max)){
  x <- as.matrix(x)
  nr <- nrow(x)
  nc <- ncol(x)
  mins <- matrix(mins, nr, nc, byrow = TRUE)
  maxs <- matrix(maxs, nr, nc, byrow = TRUE)
  out <- 2*(x - mins)/(maxs - mins) - 1
  return(out)
}
#' @rdname stdrange
ustdrange <- function(x, mins, maxs){
  x <- as.matrix(x)
  nr <- nrow(x)
  nc <- ncol(x)
  mins <- matrix(mins, nr, nc, byrow = TRUE)
  maxs <- matrix(maxs, nr, nc, byrow = TRUE)
  out <- mins + 0.5*(1 + x) * (maxs - mins)
  return(out)
}