#' Plot Track
#' 
#' Default method for plotting tracks
#' @param x object of class "track" to plot
#' @param pch default point style (filled circle)
#' @param col default color (transluscent grey)
#' @param ... other arguments to pass to \code{plot}
#' @method plot track
#' @examples
#' data(Simp)
#' is(Simp)
#' plot(Simp)

plot.track <- function(x, pch=19, col =rgb(0,0,0,.2), ...)
{
  track<-x
  x <- track$X
  y <- track$Y
  plot(x, y, asp=1, type="o", pch=pch, cex=0.5, col=col, ...)
  points(x[1], y[1], bg="green", pch=21)
  points(x[length(x)], y[length(x)], bg="red", pch=23)  
}
