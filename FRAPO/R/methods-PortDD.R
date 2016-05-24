##
## Methods for objects of class PortDD
##
## DrawDowns-method
##
setMethod(f = "DrawDowns", signature(object = "PortDD"), definition = function(object){
  ans <- slot(object, "DrawDown")
  return(ans)
})
##
## plot-method
##
setMethod(f = "plot", signature(x = "PortDD", y = "missing"), definition = function(x, main = NULL, xlab = NULL, ylab = NULL, col = c("black", "red"), grid = TRUE, invert = TRUE, ...){
  dd <- DrawDowns(x) * 100
  if(is.null(main)){main <- "Portfolio Draw Downs"}
  if(is.null(xlab)){xlab <- ""}
  if(is.null(ylab)){ylab <- "Draw Downs (percentages)"}
  if(length(col) < 2) col <- c(col[1], col[1])
  if(class(x) == "PortMdd") level <- x@MaxDD
  if(class(x) == "PortAdd") level <- x@AveDD
  if(class(x) == "PortCdd") level <- x@CDaR
  if(invert){
    dd <- -1.0 * dd
    level <- -1.0 * level
  }
  plot(dd, main = main, xlab = xlab, ylab = ylab, col = col[1], ...)
  if(grid) grid()
  abline(h = 0, col = col[1])
  abline(h = level * 100, col = col[2])
})
