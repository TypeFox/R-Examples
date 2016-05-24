#' plot_add adds a plot to a previous one
#' @title Add a plot to a previous one
#' @author Marc Girondot
#' @return Nothing
#' @param ... Parameters for plot()
#' @family plot and barplot functions
#' @description To plot data, just add use it as a normal plot. It will plot 
#' the new data without axes, or labels for axes.\cr
#' This function is complementary to matlines() and matpoints() from package graphics.
#' @examples
#' plot(x=1:100, y=sin(1:100), type="l", bty="n", xlim=c(1,200), xlab="x", ylab="y")
#' plot_add(x=1:200, y=cos(1:200), type="l", bty="n", col="red")
#' @export


plot_add <- function(...) 
  {
  	par(new=TRUE)
	par.plot <- list(...)

	par.plot <- modifyList(par.plot, list(xlab="", ylab="", main="", axes=FALSE, 
	xlim= ScalePreviousPlot()$xlim[1:2], ylim=ScalePreviousPlot()$ylim[1:2]))
	do.call(plot, par.plot) 
  }
