"identify.escouf" <-
function(x, lhorz=TRUE, lvert=TRUE, lvars=TRUE, col=2, lty=2, ...) {
	# We suppose both the RV graph and the RV.diff graph are drawn
	# So, we will use points of both graphs!
	# Calculate RV'
	n <- length(x$RV) - 1
	RVd <- x$RV[2:(n+1)] - x$RV[1:n]
	# Scale RV' to the same range as RV
	RVds <- (RVd-min(RVd))/max(RVd)*(max(x$RV)-min(x$RV))+min(x$RV)
	Xcoords <- (1:n)+0.5
	pos <- identify(c(Xcoords, Xcoords), c(x$RV[1:n], RVds), n=1, plot=FALSE)
	if (length(pos)==0)	# Operation aborted!
		stop("No position indicated on the graph!")
	if (pos>n) pos <- pos-n		# We clicked on the RV.diff curve
	# Now we calculate the corresponding level
	level <- (x$RV[pos]+x$RV[pos+1])/2
	cat("Level:", level, "\n")
	# And eventually draw the lines
	lines.escouf(x, level, lhorz, lvert, lvars, col, lty, ...)
	# We return the level obtained
	invisible(level)
}
