#' Plots object from logrr of class \code{logrrenv}. 
#'
#' @param x An object of class logrrenv to be plotted.
#' @param ... Additional graphical parameters passed to the \code{spatstat::image.im} function.  See details.
#' @param conlist Additional argument passed to the \code{spatstat::contour.im} function.
#' @param main A main title for the plot.  Default is blank.
#' @details An important aspect of this plot is the color argument (\code{col}) used for displaying the regions outside the tolerance envelopes.  If NULL (the implicit default), then the default color palette used by \code{spatstat::image.im} will be used.  Simpler schemes, e.g., c("blue", "white", "orange") can suffice. See the examples.
#' @importFrom graphics plot
#' @method plot logrrenv
#' @seealso \code{\link[spatstat]{plot.im}}
#' @seealso \code{\link[spatstat]{contour.im}}
#' @export
#' @examples
#' data(grave)
#' rsim = logrr(grave, nsim = 9)
#' plot(rsim)
#' # no border or ribben (legend).  Simple color scheme.
#' plot(rsim, col = c("blue", "white", "orange"), ribbon = FALSE, box = FALSE) 
#' # alternate color scheme
#' plot(rsim, col = topo.colors(12))
plot.logrrenv = function(x, ..., conlist = list(), main = "")
{
  if(!is.element("logrrenv", class(x))) stop("x should be an object from logrr function")
  
  # if there were no simulations
  if(is.null(x$tolenv))
  {
    spatstat::image.im(x, ...)
  }else
  {
    # create temporary im object for plotting
    xtemp = spatstat::im(mat = x$v, xcol = x$xcol, yrow = x$yrow)
    # determine which loccations within tolerance envelopes or are NA
    which_na = rbind(which(x$tolenv$v == 0, arr.ind = TRUE),
                     which(is.na(x$v), arr.ind = TRUE))

    # arguments for contour.im function
    argc = list(x = xtemp, add = TRUE, conlist)
    
    # NA any locations outisde window or inside tolerance intervals
    xtemp$v[which_na] = NA

    # arguments for image.im function
    argi = list(x = xtemp, main = main, ...)
       
    # plot colors for regions outside tolerance envelopes
    do.call(spatstat::image.im, argi)
    # add contour to plot
    do.call(spatstat::contour.im, argc) 
  }
}