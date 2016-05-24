#' Quickly finds the k'th largest element in a (large) vector
#' 
#' @param x Numeric vector.
#' @param k Integer with the rank of the desired element.
#' @return numeric k'th largest element
#' @examples y <- rnorm(1e3)
#'            stopifnot(identical(y[order(y,decreasing=TRUE)[100]],
#'                                find.kth.element(y,100)))
#' @export
find.kth.element <- function(x,k){
  Zstl_nth_element(x,k-1)[k]
}
NULL

