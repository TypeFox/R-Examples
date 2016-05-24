"extract.tsd" <-
function(e, n, series=NULL, components=NULL, ...) {
	if (missing(n)) n <- length(e$ts)
	ns <- length(e$ts)
	# if series is provided, we use it in priority
	if (is.null(series)) {
		if (n > ns) {
			warning(paste("Only", ns, "series exist in the object. Extract all series."))
			n <- ns
		}
		# We create a series value that correspond to the extraction of n first series
		series <- 1:n
	} else {					# If series is provided, we test it 
		if (is.character(series)) {
			names <- e$ts
			series <- pmatch(series, names, nomatch=0)
		} else {
			if (sum(series) > 0) series <- match(series, 1:ns, nomatch=0) else series <- c(1:ns)[match(-1:-ns, series, nomatch=0) == 0]
		}
		series <- series[series != 0]
		if (length(series) < 1)
			stop("series argument is invalid, or series does not exist in this object")
	}
	# Extract the series
	if (length(series) == 1) {
		if (ns == 1) {
			if (is.null(components)) res <- e$series else {
				if (is.character(components)) {
					names <- dimnames(e$series)[[2]]
					comp <- pmatch(components, names, nomatch=0)
				} else {
					if (sum(components) > 0) comp <- match(components, 1:ncol(e$series), nomatch=0) else comp <- c(1:ncol(e$series))[match(-1:-ncol(e$series), components, nomatch=0) == 0]
				}
				comp <- comp[comp != 0]
				if (length(comp) < 1)
					stop("No such components in the series")
				res <- e$series[, comp]
			}
		} else {
			if (is.null(components)) res <- e$series[[series]] else {
				if (is.character(components)) {
					names <- dimnames(e$series[[series]])[[2]]
					comp <- pmatch(components, names, nomatch=0)
				} else {
					if (sum(components) > 0) comp <- match(components, 1:ncol(e$series[[series]]), nomatch=0) else comp <- c(1:ncol(e$series[[series]]))[match(-1:-ncol(e$series[[series]]), components, nomatch=0) == 0]
				}
				comp <- comp[comp != 0]
				if (length(comp) < 1)
					stop("No such components in the series")
				res <- e$series[[series]][, comp]
			}
		}
	} else {
		res <- NULL
		for (i in series) {
			if (is.null(components)) ser <- e$series[[i]] else {
				if (is.character(components)) {
					names <- dimnames(e$series[[i]])[[2]]
					comp <- pmatch(components, names, nomatch=0)
				} else {
					if (sum(components) > 0) comp <- match(components, 1:ncol(e$series[[i]]), nomatch=0) else comp <- c(1:ncol(e$series[[i]]))[match(-1:-ncol(e$series[[i]]), components, nomatch=0) == 0]
				}
				comp <- comp[comp != 0]
				if (length(comp) > 0) {
					ser <- e$series[[i]][, comp]
					names <- dimnames(e$series[[i]])[[2]]
					names <- names[comp]
					if (is.null(res)) {
						res <- ser
						cnames <- paste(e$ts[i], ".", names, sep="")
					} else {
						res <- cbind(res, ser)
						cnames <- c(cnames, paste(e$ts[i], ".", names, sep=""))
					}
				}
			}
		}
		if (is.null(res))
			stop("nothing to extract!")
		dimnames(res)[[2]] <- cnames
	}
	res <- as.ts(res)
	attr(res, "units") <- e$units
	res
}
