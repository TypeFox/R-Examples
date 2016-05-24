#' Compute the bounding window from coordinates
#' 
#' @param x Columnwise matrix of coordinates, n-dimensional.
#' 
#' @details
#' Rectangular or cuboidal bounding box around coordinates given as
#' matrix 'x'. Enlarged with the Ripley and Rasson 1977 method.
#' 
#' @export

bounding_box_xy <- function(x) {
  n <- nrow(x)
  f <- (n+1)/(n-1)
  bb0 <- apply(x, 2, range)
  ce <- apply(bb0, 2, mean)
  bb1 <- t(t(bb0)-ce)
  bb2  <- bb1 * f
  t(t(bb2)+ce)
}

