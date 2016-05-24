#' pairwise distances
#' 
#' @param x Point pattern data in some format.
#' @param toroidal For toroidal distances.
#'
#' @return
#' Value is the lower triangle of the distasnce matrix, returned as a vector.
#' 
#' @useDynLib SGCS
#' @export

pairwise_distances <- function(x, toroidal=FALSE) {
  x <- internalise_pp(x)
  .External("pairwise_distances_c", x, as.integer(toroidal), package="SGCS")
}