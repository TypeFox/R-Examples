"tseries" <-
function(x) {
	if (is.null(class(x)) && class(x) != "regul" && class(x) != "tsd")
		stop("x must be a 'regul' or a 'tsd' object")
	if (class(x) == "regul") {
		if (ncol(x$y) == 1) y <- x$y[[1]] else y <- as.matrix(x$y)
		# We create a 'ts' object
		res <- ts(y, start=x$tspar$start, frequency=x$tspar$frequency)
		attr(res, "units") <- x$units
	}
	if (class(x) == "tsd") {
		if (length(x$ts) == 1) {		# We have a decomposition of a single series
			# x$series is already a ts or rts object
			res <- x$series
		} else {						# We have the decomposition of several series
			# bind all series together
			res <- x$series[[1]]
			cnames <- paste(x$ts[1], ".", dimnames(x$series[[1]])[[2]], sep="")
			for (i in 2:length(x$ts)) {
				res2 <- x$series[[i]]
				cnames <- c(cnames, paste(x$ts[i], ".", dimnames(x$series[[i]])[[2]], sep=""))
				res <- cbind(res, res2)
			}
			dimnames(res)[[2]] <- cnames
		}
		attr(res, "units") <- x$units
	}
	res
}
