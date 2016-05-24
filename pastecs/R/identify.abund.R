"identify.abund" <-
function(x, label.pts=FALSE, lvert=TRUE, lvars=TRUE, col=2, lty=2, ...) {
	p <- length(x$p.log.ind)
	Xcoords <- (1:p)
	if (label.pts == FALSE) {			# We want to identify a break
		n <- identify(c(Xcoords, Xcoords, Xcoords), c(x$p.log.ind, x$cumsum, x$p.nonull), n=1, plot=FALSE)
		if (length(n)==0)	# Operation aborted!
			stop("No position indicated on the graph!")
		if (n > p) n <- n - p		# We didn't clicked on the first curve
		if (n > p) n <- n - p		# We didn't clicked on the second curve either
		cat("Number of variables extracted:", n, "on a total of", length(x$p.log.ind), "\n")
		# And eventually draw the lines
		lines.abund(x, n, lvert, lvars, col, lty, ...)
	} else {						# We just want to label points in the graph
		vr <- as.character(x$vr)
		labels <- c(vr, vr, vr)
		n <- identify(c(Xcoords, Xcoords, Xcoords), c(x$p.log.ind, x$cumsum, x$p.nonull), labels=labels)
		n[n > p ] <- n[n > p] - p			# We didn't clicked on the first curve
		n[n > p ] <- n[n > p] - p			# We didn't clicked on the second curve either
		n <- x$vr[unique(n)]
		print(n)
	}
	# We return the value n obtained
	invisible(n)
}
