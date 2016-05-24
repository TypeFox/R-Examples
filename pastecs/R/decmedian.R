"decmedian" <-
function(x, type="additive", order=1, times=1, ends="fill") {
	call <- match.call()
	x <- as.ts(x)
	if (is.matrix(x) && ncol(x) != 1) 
		stop("only univariate series are allowed")
	if (!is.numeric(order) || order <= 0)
		stop("order must be a positive number")
	if (!is.numeric(times) || times <= 0)
		stop("times must be a positive number")
	# Check the type argument
	TYPES <- c("additive", "multiplicative")
		typeindex <- pmatch(type, TYPES)
		if (is.na(typeindex)) 
			stop("invalid type value")
		if (typeindex == -1) 
			stop("ambiguous type value")
		# make sure type is fully spelled
		type <- switch(typeindex,
				"additive"="additive",
				"multiplicative"="multiplicative")
	# Check the ends argument and treat the series accordingly
	ENDS <- c("NAs", "fill")
	endsindex <- pmatch(ends, ENDS)
	if (is.na(endsindex)) 
		stop("invalid ends value")
	if (endsindex == -1) 
		stop("ambiguous ends value")
	# make sure ends is fully spelled
	ends <- switch(endsindex,
				"NAs"="NAs",
				"fill"="fill")
	if (endsindex == 1) na.rm <- FALSE else na.rm <- TRUE
	# create our own specs component
	specs <- list(method="median", type=type, order=order, times=times, ends=ends)
	# we recuperate units from x
	units <- attr(x, "units")
	# perform filtering
	filtmedian <- function(x, n, order, term, na.rm) {
		X <- NULL
		X[(1:n) + order] <- x
		X[1:order] <- NA
		X[(1:order)+n+order] <- NA
		f <- NULL
		for (i in (1:n))
			f[i] <- median(X[(1:term) + i - 1], na.rm=na.rm)
		f
	}
	term <- 2*order + 1
	n <- length(x)
	filtered <- x						# We don't change the initial series, but a copy of it
	for (i in 1:times)
		filtered <- filtmedian(filtered, n=n, order=order, term=term, na.rm=na.rm)
	filtered <- ts(filtered, start=start(x), frequency=frequency(x))
	# Calculate residuals
	if (type == "additive") {
		residuals <- x - filtered
	} else {
		residuals <- x / filtered
	}
	series <- ts.union(filtered, residuals)
	# create our own 'tsd' structure
	res <- list(ts="series", series=series, units=units, specs=specs, call=call)
	class(res) <- "tsd"		# change the class of the object to 'tsd'
	res
}
