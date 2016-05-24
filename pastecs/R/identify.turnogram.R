"identify.turnogram" <-
function(x, lvert=TRUE, col=2, lty=2, ...) {
	n <- length(x$interval)
	if (x$type == "Complete") {
		pos <- identify(c(x$interval, x$interval, x$interval), c(x$info, x$info.min, x$info.max), n=1, plot=FALSE)
	} else {
		pos <- identify(x$interval, x$info, n=1, plot=FALSE)
	}
	if (length(pos)==0)	# Operation aborted!
		stop("No position indicated on the graph!")
	if (pos > n) pos <- pos - n		# We clicked on either min or max curve
	if (pos > n) pos <- pos - n		# We clicked on max curve
	# Now we calculate the corresponding level
	level <- x$interval[pos]
	cat("Level      :", level, "\n")
	cat("Information:", x$info[pos], "\n")
	cat("Probability:", 2^-abs(x$info[pos]), "\n")
	cat("Nbr of obs.:", x$n[pos], "\n")
	cat("Turnpoints :", x$turns[pos], "\n")
	# And eventually draw the lines
	if (lvert == TRUE) abline(v=level, col=col, lty=lty, ...)
	# We return the level obtained
	invisible(level)
}
