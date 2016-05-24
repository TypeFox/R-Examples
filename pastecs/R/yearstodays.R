"yearstodays" <-
function (x, xmin=NULL) {
	x <- x
	xmin <- xmin
	if (!is.null(xmin)) {		# xmin is given
		x <- x * 365.25
		x <- x - min(x, na.rm=TRUE) + xmin
	} else {					# xmin is not given, we construct julian dates (compatibles with chron, dates,...)
		if(is.null(yearorig <- options("chron.origin")$year))
			yearorig <- 1970
		x <- (x - yearorig) * 365.25
	}
	x
}
