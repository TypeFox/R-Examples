#' Extract clusters from clustering object.
#'
#' @param x object
#' @export
#' @keywords internal
clusters <- function(x) UseMethod("clusters", x)
#' @export
clusters.kmeans <- function(x) as.vector(x$cluster)
#' @export
clusters.default <- function(x) as.vector(x)
#' @export
clusters.partition <- function(x) x$clustering
