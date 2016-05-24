#' Impute Optimal Scores for Rating Categories
#' 
#' Replace original ratings with optimal scores based on \code{\link{cds}} output..
#' 
#' @param object An object of class \code{cds}
#' @param data An object of class \code{cdsdata} to be cleaned, or the original data.
#' @param K The number of classes in the solution that must be kept.
#' @param col.subset An optional subset
#' @param \dots Additional arguments.
#' @keywords utility
#' @export 
clean.scales <- function(object, data, K, col.subset = NULL, ...) {
  UseMethod("clean.scales")  
}
