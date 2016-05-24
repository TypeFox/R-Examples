#' ggtern Constructor
#' 
#' Plots in \code{ggtern} are instigated via the default constructor: \code{ggtern(...)}, which is essentially a convenience wrapper for the following: 
#' \code{ggplot{...} + coord_tern()}, indeed, if one wishes to use \code{ggplot{...} + coord_tern()} then this is quite satisfactory.
#' @inheritParams ggplot2::ggplot
#' @aliases constructor
#' @return \code{ggtern(...)} returns an object of class \code{ggplot}.
#' @seealso For an introduction to the \code{ggtern} package, (including many examples), click \link[=ggtern-package]{HERE}.
#' @export
ggtern <- function(...){
  ggplot(...) + coord_tern()
}