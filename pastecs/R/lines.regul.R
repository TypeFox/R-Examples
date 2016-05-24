"lines.regul" <-
function(x, series=1, col=3, lty=1, plot.pts=TRUE,...) {
	# The next function actually draw the lines
	regul.lines <- function(X, Series, Col, Lty, Plot.pts, ...) {
		i <- Series
		if (i > ncol(X$y))
			stop("This series does not exist")
		# Calculate the time span
		xlbr <- min(X$x, na.rm=TRUE)
		xubr <- max(X$x, na.rm=TRUE)
		# Calculate the y span
		ylbr <- min(X$y[,i], na.rm=TRUE)
		yubr <- max(X$y[,i], na.rm=TRUE)
		# Trace the regulated series (but without NA values)
		xv <- X$x
		yv <- X$y[,i]
		xv <- xv[!is.na(yv)]
		yv <- yv[!is.na(yv)]
		lines(xv, yv, col=Col, lty=Lty)
		if (Plot.pts == TRUE) {					# plot points of regular series
			points(xv, yv, col=Col, pch="+")
			# Indicate matching points
			points(X$x[is.finite(X$match.dist)], X$y[is.finite(X$match.dist), i], col=Col, pch="O")
		}
		# Indicate spanning of regulated series
		lines(c(xlbr, xlbr), c(ylbr, yubr/3*2), col=Col, lty=2, type="l")
		lines(c(xubr, xubr), c(ylbr, yubr/3*2), col=Col, lty=2, type="l")
	}
	invisible(regul.lines(x, series, col, lty, plot.pts, ...))
}
