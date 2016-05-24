# plot method for cancor, borrowing from interpPlot

default.arg <- function(args.list, arg, default){
    if (is.null(args.list[[arg]])) default else args.list[[arg]]
}

plot.cancor <- function(x, which=1, 
		xlim, ylim, xlab, ylab,
    points=TRUE, add=FALSE, col=palette()[1],
		ellipse = TRUE, ellipse.args = list(), 
    smooth=FALSE, smoother.args = list(), col.smooth=palette()[3],
		abline=TRUE, col.lines = palette()[2], lwd=2,
		labels=rownames(xy), 
		id.method = "mahal", 
		id.n = 0, id.cex = 1, id.col = palette()[1],
		...) {

  if (!inherits(x, "cancor")) 
      stop("Not a cancor object")
  # sanity checks on which
  if (length(which)>1) {
  	which <- which[1]
  	warning("plot.cancor only plots one dimension at a time")
  }
  if (! which %in% 1:x$ndim) stop(paste(which, "is not among the canonical dimensions"))

	xy <- cbind(scores(x, type="x")[,which],
              scores(x, type="y")[,which])
	
	lims <- apply(xy, 2, range)
	if (missing(xlim)) xlim=lims[,1]
	if (missing(ylim)) ylim=lims[,2]
  
  if (missing(xlab)) xlab <- paste(x$names$set.names[1], "dimension", which) 
  if (missing(ylab)) ylab <- paste(x$names$set.names[2], "dimension", which) 
    
	if (!add) plot(xy, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)  
	if (points) points(xy, col=col, ...)
	
  if (ellipse) {
			levels <- default.arg(ellipse.args, "levels", 0.68)
      ellipse.args <- c(list(xy, add = TRUE, plot.points = FALSE, levels=levels), 
          ellipse.args)
      do.call(dataEllipse, ellipse.args)
  }
  if (abline) {
    abline(lsfit(xy[, 1], xy[, 2]), col = col.lines, 
        lwd = lwd)
  }
  if (!is.null(labels)) {
  showLabels(xy[, 1], xy[, 2], labels = labels, id.method = id.method, 
        id.n = id.n, id.cex = id.cex, id.col = id.col)
  }
  if (smooth) {
    loessLine(xy[,1], xy[,2], col=col.smooth, smoother.args = smoother.args, log.x=FALSE, log.y=FALSE)
  }
#  invisible()
}
