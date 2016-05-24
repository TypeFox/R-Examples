"decdiff" <-
function(x, type="additive", lag=1, order=1, ends="fill") {
	call <- match.call()
	x <- as.ts(x)
	if (is.matrix(x) && ncol(x) != 1) 
		stop("only univariate series are allowed")
	if (!is.numeric(lag) || lag <= 0)
		stop("lag must be a positive number")
	if (!is.numeric(order) || order <= 0)
		stop("order must be a positive number")
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
	ENDS <- c("NAs", "fill", "drop")
	endsindex <- pmatch(ends, ENDS)
	if (is.na(endsindex)) 
		stop("invalid ends value")
	if (endsindex == -1) 
		stop("ambiguous ends value")
	# make sure ends is fully spelled
	ends <- switch(endsindex,
				"NAs"="NAs",
				"fill"="fill",
				"drop"="drop")
	# create our own specs component
	specs <- list(method="diff", type=type, lag=lag, order=order, ends=ends)
	# we recuperate units from x
	units <- attr(x, "units")
	# The next function add enough data to the left (either NA or the mean of first few values)
	# to obtain a series of the same length as x after difference
	padleft <- function(x, Lag, fill) {
		x <- window(x, start=start(lag(x, Lag)), end=end(x), extend=TRUE)
		if (fill == TRUE)			# We fill padded data with the mean of first few values
			x[1:Lag] <- mean(x[(1:Lag)+Lag], na.rm=TRUE)
		x
	}
	filtered <- switch(endsindex,
					"NAs"=padleft(x, lag*order, fill=FALSE),				# We add NA's in front of the series
					"fill"=padleft(x, lag*order, fill=TRUE),				# We add the mean of first values in front of the series
					"drop"=x)												# We keep x like that
	# perform filtering
	filtered <- diff(filtered, lag=lag, difference=order)
	# Calculate residuals
	if (type == "additive") {
		residuals <- x - filtered
	} else {
		residuals <- x / filtered
	}
	series <- ts.intersect(filtered, residuals)
	# create our own 'tsd' structure
	res <- list(ts="series", series=series, units=units, specs=specs, call=call)
	class(res) <- "tsd"		# change the class of the object to 'tsd'
	res
}
