#' @name plotn
#' @author Sven E. Templer
#' @title Plot Nothing (but a Plot)
#' @description 
#' Create a plot, with empty elements by presetting default parameters.
#' @details
#' For details about the function see \link{plot}, which is called from
#' \code{plotn}. More detailed information in \link{plot.default} and
#' \link{par}.
#' @param x Coordinates of the points.
#' @param y Coordinates of the y-axis.
#' @param type Plot type.
#' @param xlab,ylab Axis titles.
#' @param xaxt,yaxt Axis types.
#' @param frame.plot Plot the frame.
#' @param ... Forwarded arguments to \code{plot}.

#' @export
plotn <- function(x = 0:1, y = NULL, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame.plot=F, ...) {
  plot(x,y,type=type,xlab=xlab,ylab=ylab,xaxt=xaxt,yaxt=yaxt,frame.plot=frame.plot,...)
}
