
##==============================================================================
##  emptyplot: open plot region for shape plotting
##==============================================================================

emptyplot <- function (xlim=c(0,1), ylim=xlim, asp=1,
  frame.plot = FALSE, col=NULL, ...) {

  plot(0, type = "n", xlab="", ylab = "", asp=asp, axes=FALSE,
       frame.plot = frame.plot,
       xlim = xlim, ylim = ylim, xaxs="i", yaxs="i", ...)

  if ( ! is.null(col) )
    rect(xlim[1], ylim[1], xlim[2], ylim[2], col=col)

}
