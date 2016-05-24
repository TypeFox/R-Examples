#' ScalePreviousPlot returns the scale of the previous plot
#' @title Return the scale of the previous plot
#' @author Marc Girondot
#' @return A list with xlim and ylim
#' @family plot and barplot functions
#' @description Return a list with the limits of the previous plot, the center, the range, and the position of label on this axe. 
#' @examples
#' par(xaxs="i", yaxs="i")
#' plot(x=1:100, y=sin(1:100), type="l", bty="n", xlim=c(1,200), xlab="x", ylab="y")
#' xlim= ScalePreviousPlot()$xlim[1:2]
#' ylim= ScalePreviousPlot()$ylim[1:2]
#' par(xaxs="r", yaxs="i")
#' plot(x=1:100, y=sin(1:100), type="l", bty="n", xlim=c(1,200), xlab="x", ylab="y")
#' xlim= ScalePreviousPlot()$xlim[1:2]
#' ylim= ScalePreviousPlot()$ylim[1:2]
#' # Here is an example of the use of the label output
#' plot(x=1:100, y=sin(1:100), type="l", bty="n", xlim=c(1,200), xlab="", ylab="")
#' text(x=ScalePreviousPlot()$xlim["label"], y=ScalePreviousPlot()$ylim["center"], 
#'   xpd=TRUE, "Legend for Y axes", pos=3, srt=90)
#' text(x=ScalePreviousPlot()$xlim["center"], y=ScalePreviousPlot()$ylim["label"], 
#'   xpd=TRUE, "Legend for X axes", pos=1)
#' @export


ScalePreviousPlot <- function() {
  if (par("xaxs")=="i") {
    x1 <- par("usr")[1]
    x2 <- par("usr")[2]
  } else {
    x2 <- (par("usr")[1]+par("usr")[2]*26)/27
    x1 <- x2*26-par("usr")[2]/0.04
  }
    if (par("yaxs")=="i") {
    y1 <- par("usr")[3]
    y2 <- par("usr")[4]
  } else {
    y2 <- (par("usr")[3]+par("usr")[4]*26)/27
    y1 <- y2*26-par("usr")[4]/0.04
  }
  return(list(xlim=c(begin=x1, end=x2, center=x1+(x2-x1)/2, range=x2-x1, label=x1-(x2-x1)/2.6), 
  				ylim=c(begin=y1, end=y2, center=y1+(y2-y1)/2, range=y2-y1, label=y1-(y2-y1)/3)))
}
