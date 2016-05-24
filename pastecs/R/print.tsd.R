"print.tsd" <-
function(x, ...) {
    cat(" Call:\n")
    cat(" ")
    dput(x$call)
    ser <- x$ts
    if (length(ser) == 1) {		# We have the decomposition of a single series in the object
    	cat("\nComponents\n")
    	print(dimnames(x$series)[[2]])
    } else {					# We have the decomposition of several series in the object
    	cat("\nSeries:\n")
		print(x$ts)
		cat("\n")
		for (i in 1:length(ser)) {
			cat("\nComponents for", ser[i], "\n")
			print(dimnames(x$series[[i]])[[2]])
    	}
	}
	invisible(x)
}
