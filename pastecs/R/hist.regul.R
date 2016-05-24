"hist.regul" <-
function(x, nclass=30, col=c(4, 5, 2), xlab=paste("Time distance in", x$units, "with start =", min(x$x), ", n = ", length(x$x), ", deltat =", x$tspar$deltat), ylab=paste("Frequency, tol =", x$specs$tol), main="Number of matching observations", plotit=TRUE, ...) {
	# The next function actually draw the histogram
	regul.hist <- function(X, Col, Xlab, Ylab, Main, PlotIt, ...) {
		# Prepare the vector of data
		if (X$specs$tol.type == "none")
			stop("tol.type was 'none', all observations were interpolated!")
		Tol <- X$specs$tol
		if (Tol == 0) HT <- 1.001 else HT <- 101*Tol/100
		Data <- abs(X$match.dist)
		Data[is.infinite(Data)] <- HT 			# Inf are replaced by a value higher than Tol
		Data[Data == 0] <- -0.00001				# For putting exact matching values in a separate category
		# Don't draw, but get vectors of results
		res <- hist(Data, nclass=nclass, plot=FALSE)
		classes <- res$breaks[2:length(res$breaks)]
		ncl <- length(classes)
		classes[ncl] <- Inf
		counts <- res$counts
		counts <- counts[counts != 0]
		lc <- length(counts)
		counts2 <- NULL
		for (i in 1:lc) {
			counts2[i] <- sum(counts[1:i])
		}
		names(counts2) <- names(counts)
		# Create a vector for colors, so as the first and last classes are drawn in a different color
		cols <- NULL
		cols[1] <- Col[1]
		if (ncl > 2) cols[2:(ncl-2)] <- Col[2]
		cols[ncl] <- Col[3]
		# Actually draw the histogram
		if (PlotIt == TRUE)
			hist(Data, nclass=nclass, col=cols, xlab=Xlab, ylab=Ylab, main=Main)
		counts2
	}
	invisible(regul.hist(x, col, xlab, ylab, main, plotit, ...))
}
