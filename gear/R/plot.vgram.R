#' Plot vgram object
#' 
#' Plots object of class '\code{vgram}' produced by the \code{vgram} function.  The plot is based on the \code{lattice::xyplot} function.
#' 
#' @param x A \code{vgram} object produced by the \code{vgram} function.
#' @param ... Additional arguments to pass the \code{lattice::xyplot} function to change aspects of the plot.
#' @param split A logical value indicating whether, for a directional semivariogram, the directional semivariograms should be displayed in a single or split panels.  Default is FALSE, for a single panel.
#' @return NULL
#' @author Joshua French
#' @export
#' @importFrom lattice xyplot
#' @method plot vgram
#' @seealso \code{\link[lattice]{xyplot}}, \code{\link{vgram}}
#' @examples 
#' data(co)
#' v = vgram(Al ~ 1, co, ~ easting + northing)
#' plot(v)
#' v2 = vgram(Al ~ 1, co, ~ easting + northing, angle = 22.5, ndir = 4)
#' plot(v2)
#' # show how attributes can be changed using different arguments
#' # available in \code{lattice::xyplot}.
#' plot(v2, col = 2:5)
#' plot(v2, col = 2:5, pch = 1:4)
#' plot(v2, col = 2:5, pch = 1:4, lty = 2:5, type = "b")
#' plot(v2, col = 2:5, pch = 1:4, lty = 2:5, type = "b", 
#'  key=list(text=list(levels(as.factor(v2$semi$angle))), 
#'  space='right', points=list(pch=1:4, col=2:5), 
#'  lines=list(col=2:5, lty = 2:5)))
#'plot(v2, split = TRUE)

plot.vgram = function(x, ..., split = FALSE)
{
  if(class(x) != "vgram") stop("x should be a vgram object produced by the vgram function")
  if(x$ndir == 1)
  {
    lattice::xyplot(semivariance ~ distance, data = x$semi, ...)
  }else
  {
    if(split)
    {
      lattice::xyplot(semivariance ~ distance | angle, 
                      data = x$semi, ...)
    }else
    {
      lattice::xyplot(semivariance ~ distance, data = x$semi, 
                      groups = x$semi$angle, ...)
    }
  }
}
