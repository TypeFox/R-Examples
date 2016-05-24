"lines.turnpoints" <-
function(x, max=TRUE, min=TRUE, median=TRUE, col=c(4, 4, 2), lty=c(2, 2, 1), ...) {
	# The next function actually draws the graph
	turnpoints.lines <- function(X, Max, Min, Median, Col, Lty, ...) {
		x.peaks <- X$pos[X$peaks]
		y.peaks <- X$points[X$peaks]
		y.peaks.approx <- approx(x.peaks, y.peaks, X$pos, method="linear")$y
		x.pits <- X$pos[X$pits]
		y.pits <- X$points[X$pits]
		y.pits.approx <- approx(x.pits, y.pits, X$pos, method="linear")$y
		y.median <- y.pits.approx + (y.peaks.approx - y.pits.approx) / 2
		if (Max)
			lines(x.peaks, y.peaks, col=Col[1], lty=Lty[1], ...)
		if (Min)
			lines(x.pits, y.pits, col=Col[2], lty=Lty[2], ...)
		if (Median)
			lines(X$pos, y.median, col=Col[3], lty=Lty[3], ...)
	}
	col <- rep(col, length.out=3)
	lty <- rep(lty, length.out=3)
	invisible(turnpoints.lines(x, max, min, median, col, lty, ...))	
}
