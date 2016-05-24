#' Plot a density surface model.
#'
#' See \code{\link{plot.gam}}.
#'
#' @aliases plot.dsm
#'
#' @param x a \code{dsm} object
#' @param \dots other arguments passed to \code{\link{plot.gam}}.
#' @return a plot!
#' @export
#'
#' @seealso dsm plot.gam
#' @author David L. Miller
#' @importFrom graphics plot
#'
plot.dsm<-function(x,...){
  class(x) <- class(x)[class(x)!="dsm"]
  return(plot(x,...))
}
