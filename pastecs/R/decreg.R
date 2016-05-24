"decreg" <-
function(x, xreg, type="additive") {
	call <- match.call()
	if (is.matrix(x) && ncol(x) != 1) 
		stop("only univariate series are allowed")
	if (length(x) != length(xreg))
		stop("x and xreg must have same row number")	
	x <- as.ts(x)
	xreg <- as.ts(xreg)
	# Make sure "tsp" attributes are the same for both series
	attr(xreg, "tsp") <- attr(x, "tsp")
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
	# create our own specs component
	specs <- list(method="reg", type=type, xreg=xreg)
	# we recuperate units from x
	units <- attr(x, "units")
	model <- xreg
	# Calculate residuals
	if (type == "additive") {
		residuals <- x - model
	} else {
		residuals <- x / model
	}
	series <- ts.union(model, residuals)
	# create our own 'tsd' structure
	res <- list(ts="series", series=series, units=units, specs=specs, call=call)
	class(res) <- "tsd"		# change the class of the object to 'tsd'
	res
}
