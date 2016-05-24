#' Plot a fitted detection function
#'
#' This is just a simple wrapper around \code{\link{plot.ds}}. See the manual
#' page for that function for more information.
#'
#' @param x an object of class \code{dsmodel}.
#' @param ... extra arguments to be passed to \code{\link{plot.ds}}.
#' @return \code{NULL}, just produces a plot.
#' @aliases plot.dsmodel
#' @export
#' @author David L. Miller
#' @importFrom graphics plot
plot.dsmodel <- function(x,...){

  plot(x$ddf,...)

  invisible()
}
