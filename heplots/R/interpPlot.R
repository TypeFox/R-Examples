# Plot an interpolation between two related data sets,
# typically transformations of each other.
# This function is designed to be used in animations.

# DONE: incorporate xlim, ylim from xy1, xy2
# TODO: incorporate the loop over alpha values?
# DONE: allow add=TRUE to add to an existing plot
# DONE: add col= arg
# DONE: return invisibly the interpolated coordinates

interpPlot <- function(xy1, xy2, alpha,
		xlim, ylim, points=TRUE, add=FALSE, col=palette()[1],
		ellipse = FALSE, ellipse.args = NULL, 
		abline=FALSE, col.lines = palette()[2], lwd=2,
		id.method = "mahal", labels=rownames(xy1), 
		id.n = 0, id.cex = 1, id.col = palette()[1],
		segments=FALSE, segment.col="darkgray",
		 ...){

	xy1 <- as.matrix(xy1)
	xy2 <- as.matrix(xy2)
	# interpolate
	xy <- xy1 + alpha * (xy2-xy1)

	# set default plot limits to incorporate xy1, xy2
	lims <- apply(rbind(xy1, xy2), 2, range)
	if (missing(xlim)) xlim=lims[,1]
	if (missing(ylim)) ylim=lims[,2]

	if (!add) plot(xy, xlim=xlim, ylim=ylim, ...)  # can't use type="n" here
	if (points) points(xy, col=col, ...)
	
  if (ellipse) {
      ellipse.args <- c(list(xy, add = TRUE, plot.points = FALSE), 
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
  
  if(segments) {
  	segments(xy1[,1], xy1[,2], xy[,1], xy[,2], col=segment.col)
  }
  invisible(xy)
}

CIRCLETEST <-FALSE

if(CIRCLETEST) {
	# two circles
	circle <- function(radius, segments=40) {
			angles <- (0:segments)*2*pi/segments
			radius * cbind( cos(angles), sin(angles))
	}
	C1 <- circle(2)
	C2 <- circle(4)
	
	# show as separate frames
	for (a in seq(0,1,.1)) {
		interpPlot(C1, C2, alpha=a, pch=16, col=rgb(1,0,0,a), type='b', asp=1)
	}

	# show all in one frame
	for (a in seq(0,1,.1)) {
		interpPlot(C1, C2, alpha=a, pch=16, col=rgb(1,0,0,a), type='b', add=a>0)
	}
}
