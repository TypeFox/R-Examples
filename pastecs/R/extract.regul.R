"extract.regul" <-
function(e, n, series=NULL, ...) {
	if (missing(n)) n <- ncol(e$y)
	nc <- ncol(e$y)
	if (is.null(nc)) nc <- 1		# if ncol() return null, we have just a single vector
	# if series is provided, we use it in priority
	if (is.null(series)) {
		if (n > nc) {
			if (nc > 1) warning(paste("Only", nc, "series exist in the object. Extract all series."))
			n <- nc
		}
		# We create a series value that correspond to the extraction of n first series
		series <- 1:n
	}
	if (nc == 1) {
		warning("Only one series in the object. Extract it.")
		y <- e$y[[1]]
	} else {			# Use series to determine which series to extract
		y <- as.matrix(e$y)[, series]
	}		
	# We create a 'ts' object
	res <- ts(y, start = e$tspar$start, frequency = e$tspar$frequency)
	attr(res, "units") <- e$units
	res
}
