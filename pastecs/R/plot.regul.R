"plot.regul" <-
function(x, series=1, col=c(1,2), lty=c(par("lty"), par("lty")), plot.pts=TRUE, leg=FALSE, llab=c("initial", x$specs$methods[series]), lpos=c(1.5, 10), xlab=paste("Time (", x$units, ")", sep=""), ylab="Series", main=paste("Regulation of", names(x$y)[series]), ...) {
	# The next function actually draw the graph
	regul.graph <- function(X, Series, Col, Lty, Plot.pts, Leg, Llab, Lpos, Xlab, Ylab, Main, ...) {
		i <- Series
		if (i > ncol(X$y))
			stop("This series does not exist")
		# Calculate the total time span
		xlbi <- min(X$xini, na.rm=TRUE)
		xubi <- max(X$xini,na.rm=TRUE)
		xlbr <- min(X$x, na.rm=TRUE)
		xubr <- max(X$x, na.rm=TRUE)
		xlb <- min(xlbi, xlbr)
		xub <- max(xubi, xubr)
		xspan <- c(xlb, xub)
		# Calculate the y span
		ylbi <- min(X$yini[,i], na.rm=TRUE)
		yubi <- max(X$yini[,i], na.rm=TRUE)
		ylbr <- min(X$y[,i], na.rm=TRUE)
		yubr <- max(X$y[,i], na.rm=TRUE)
		ylb <- min(ylbi, ylbr)
		yub <- max(yubi, yubr)
		yspan <- c(ylb, yub)
		plot(xspan, yspan, type="n", xlab=Xlab, ylab=Ylab, main=Main, ...)
		# Trace the initial series
		lines(X$xini, X$yini[,i], col=Col[1], lty=Lty[1])
		# Trace the regulated series (but without NA values)
		xv <- X$x
		yv <- X$y[,i]
		xv <- xv[!is.na(yv)]
		yv <- yv[!is.na(yv)]
		lines(xv, yv, col=Col[2], lty=Lty[2])
		if (Plot.pts == TRUE) {					# plot points of regular series
			points(xv, yv, col=Col[2], pch="+")
			# Indicate matching points
			points(X$x[is.finite(X$match.dist)], X$y[is.finite(X$match.dist), i], col=Col[2], pch="O")
		}
		# Indicate respective spanning of initial and regulated series
		lines(c(xlbi, xlbi), c(ylb+yub/3, yub), col=Col[1], lty=2, type="l")
		lines(c(xubi, xubi), c(ylb+yub/3, yub), col=Col[1], lty=2, type="l")
		lines(c(xlbr, xlbr), c(ylb, yub/3*2), col=Col[2], lty=2, type="l")
		lines(c(xubr, xubr), c(ylb, yub/3*2), col=Col[2], lty=2, type="l")
		# If Leg is TRUE, print a legend
		if (Leg == TRUE) {
			legend(Lpos[1], Lpos[2], Llab, col=Col, lty=Lty)
		}
	}
	invisible(regul.graph(x, series, col, lty, plot.pts, leg, llab, lpos, xlab, ylab, main, ...))
}
