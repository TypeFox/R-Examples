#' Lines method for Bolstad objects
#' 
#' Allows simple addition of posterior distributions from other results to an
#' existing plot
#' 
#' @param x an object of class \code{Bolstad}.
#' @param \dots any additional parameters to be passed to \code{graphics::lines}.
#' 
#' @method lines Bolstad
#' 
#' @export

lines.Bolstad = function(x, ...)
  lines(x$param.x, x$posterior, ...)
