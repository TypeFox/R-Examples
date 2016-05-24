#' Plots object from kdest of class \code{kdenv}. 
#'
#' @param x An object of class kdenv to be plotted.
#' @param ... Additional graphical parameters passed to \code{plot.fv} function.
#' @param shadecol1 Shade color for max/min envelopes. 
#' @param shadecol2 Shade color for confidence envelopes.
#' @param main A main title for the plot.  Default is blank.
#' @param legend Logical for whether a legend should automatically be displayed.  Default if \code{FALSE}.
#' @import spatstat
#' @importFrom graphics plot polygon
#' @method plot kdenv
#' @seealso \code{\link[spatstat]{plot.fv}}
#' @export
#' @examples
#' data(grave)
#' kd1 = kdest(grave, nsim = 19, level = 0.9)
#' plot(kd1)

plot.kdenv = function(x, ..., shadecol1 = "grey", shadecol2 = "lightblue", main = "", legend = FALSE)
{
  if(!is.element("kdenv", class(x))) stop("x should be an object from kdest function")
  # if there were no simulations
  if(length(x) == 1) spatstat::plot.fv(x[[1]], main = main, legend = legend, ...)
  if(length(x) > 1)
  {
    # create main plot
    spatstat::plot.fv(x[[1]], main = main, legend = legend, shadecol = shadecol1, ...)
    # shade confidence bands
    # do some additional prep in case xlim is specified, otherwise
    # the polygon will go beyond the desired boundary
    plotargs = list(...)
    xmax = max(x$r)
    if(!is.null(plotargs$xlim)){ xmax = max(plotargs$xlim) }
    pin = which(x$r <= xmax)
    
    graphics::polygon(c(x$r[pin], rev(x$r[pin])), c(x$qhi[pin], rev(x$qlo[pin])), 
                      col = shadecol2, border = NA)
    # redraw parts of plot covered by confidence bands
    spatstat::plot.fv(x[[1]], main = "", legend = FALSE, shadecol = NA, add = TRUE, ...)
  }
}