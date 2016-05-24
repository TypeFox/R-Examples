#' barplot_errbar plot a barplot with error bar on y
#' @title Plot a barplot graph with error bar on y
#' @author Marc Girondot
#' @return A numeric vector (or matrix, when beside = TRUE), say mp, giving the coordinates of all the bar midpoints drawn, useful for adding to the graph.\cr
#' If beside is true, use colMeans(mp) for the midpoints of each group of bars, see example.
#' @param ... Parameters for barplot() such as main= or ylim=
#' @param errbar.y The length of error bars for y. Recycled if necessary.
#' @param errbar.y.plus The length of positive error bars for y. Recycled if necessary.
#' @param errbar.y.minus The length of negative error bars for y. Recycled if necessary.
#' @param y.plus The absolut position of the positive error bar for y. Recycled if necessary.
#' @param y.minus The absolut position of the nagative error bar for y. Recycled if necessary.
#' @param errbar.tick Size of small ticks at the end of error bars defined as a proportion of total width or height graph size.
#' @param errbar.lwd Error bar line width, see par("lwd")
#' @param errbar.lty Error bar line type, see par("lwd")
#' @param errbar.col Error bar line color, see par("col")
#' @param add If true, add the graph to the previous one.
#' @family plot and barplot functions
#' @description To plot data, just use it as a normal barplot but add the errbar.y 
#' values or errbar.y.minus, errbar.y.plus if bars for y axis are 
#' asymetric. Use y.plus and y.minus to set absolut limits for
#' error bars. Note that y.plus and y.minus have priority over errbar.y, 
#' errbar.y.minus and errbar.y.plus.
#' @seealso \code{plot_errorbar}
#' @examples
#' barplot_errbar(rnorm(10, 10, 3), 
#'		xlab="axe x", ylab="axe y", bty="n", 
#' 		errbar.y.plus=rnorm(10, 1, 0.1), col=rainbow(10), 
#' 		names.arg=paste("Group",1:10), cex.names=0.6)
#' y <- rnorm(10, 10, 3)
#' barplot_errbar(y, 
#'                	xlab="axe x", ylab="axe y", bty="n", 
#'             		y.plus=y+2)
#' @export


barplot_errbar <- function(..., 
                        errbar.y=NULL, 
                        errbar.y.plus=NULL, errbar.y.minus=NULL,
                        y.plus=NULL, y.minus=NULL,
                        errbar.tick=1/50, 
                        errbar.lwd=par("lwd"), 
                        errbar.lty=par("lty"), 
                        errbar.col=par("fg"), 
                        add=FALSE) 
  {
# errbar.y=NULL; errbar.y.plus=NULL; errbar.y.minus=NULL; y.plus=NULL; y.minus=NULL; errbar.tick=1/50; errbar.lwd=par("lwd"); errbar.lty=par("lty"); errbar.col=par("fg"); add=FALSE 

  par.plot <- list(...)
  if (add) {
  	s <- ScalePreviousPlot()
  	par(new=TRUE)
  	par.plot <- modifyList(par.plot, list(xlim=s$xlim, ylim=s$ylim, xlab="", ylab="", main="", axes=FALSE))
  }


  essai <- do.call(barplot, par.plot)
  x <- as.numeric(essai) 
  y <- par.plot[[1]]
  
  if (!is.null(y.plus)) errbar.y.plus <- y.plus-y
  if (!is.null(y.minus)) errbar.y.minus <- y-y.minus
  
  if (is.null(errbar.y.minus) & !is.null(errbar.y)) {
  	errbar.y.minus <- errbar.y
  }
  if (is.null(errbar.y.plus) & !is.null(errbar.y)) {
  	errbar.y.plus <- errbar.y
  }
    
  sizebar <- (par("usr")[2]-par("usr")[1])*errbar.tick
  
  par(xpd=TRUE)
  
  if (!is.null(errbar.y.minus)) {
    segments(x, y-errbar.y.minus, x, y, 
             col=errbar.col, lty=errbar.lty, lwd=errbar.lwd)
    segments(x-sizebar, y-errbar.y.minus, x+sizebar, y-errbar.y.minus, 
             col=errbar.col, lty=errbar.lty, lwd=errbar.lwd)
  }
  if (!is.null(errbar.y.plus)) {
    segments(x, y+errbar.y.plus, x, y, 
             col=errbar.col, lty=errbar.lty, lwd=errbar.lwd)
    segments(x-sizebar, y+errbar.y.plus, x+sizebar, y+errbar.y.plus, 
             col=errbar.col, lty=errbar.lty, lwd=errbar.lwd)
  }

  return(invisible(essai))
}
