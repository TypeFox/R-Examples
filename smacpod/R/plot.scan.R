#' Plots object from spscan.test of class \code{scan}. 
#'
#' @param x An object of class scan to be plotted.
#' @param ... Additional graphical parameters passed to the \code{spatstat::plot.ppp} function.
#' @param nv The number of verticies to draw the circle. Default is 100.
#' @param border The border color of the circle.  Default is NULL, meaning black.
#' @param ccol Fill color of the circles.  Default is NULL, indicating empty.
#' @param clty Line type of circles.  Default is NULL, indicting lty = 1.
#' @param clwd Line width of circles.  Default is NULL, indicating lwd = 2 for the most likely cluster and lwd = 1 for the rest.  
#' @importFrom graphics plot
#' @importFrom plotrix draw.circle
#' @method plot scan
#' @seealso \code{\link[spatstat]{plot.ppp}}
#' @export
#' @examples
#' data(grave)
#' out = spscan.test(grave, case = 2, alpha = 0.1)
#' plot(out, chars = c(1, 20), main = "most likely cluster")
plot.scan = function(x, ..., nv = 100, border = NULL, ccol = NULL, clty = NULL, clwd = NULL)
{
  if(class(x) != "scan") stop("x should be a scan object from spscan.test function")
  spatstat::plot.ppp(x$ppp, ...)
  
  # number of centroids
  nc = length(x$clusters)
  
  # set default values
  if(is.null(border)) border = rep(1, nc)
  if(is.null(clty)) clty = rep(1, nc)
  if(is.null(ccol)) ccol = rep(NA, nc)
  if(is.null(clwd)) clwd = c(2, rep(1, nc - 1))
  
  # more sanity checking
  if(length(border) != nc) stop("if specific, border must have length equal to nrow(x$coords)")
  if(length(clty) != nc) stop("if specific, border must have length equal to nrow(x$coords)")
  if(length(ccol) != nc) stop("if specific, border must have length equal to nrow(x$coords)")
  if(length(clwd) != nc) stop("if specific, border must have length equal to nrow(x$coords)")
  
  # plot clusters
  for(i in 1:nc)
  {
    plotrix::draw.circle(x$clusters[[i]]$coords[1, 1], x$clusters[[i]]$coords[1, 2], 
                         x$clusters[[i]]$r, nv = nv, border = border, col = ccol[i], 
                         lty = clty[i], lwd = clwd[i])
  }
}
