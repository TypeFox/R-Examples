"decloess" <-
function(x, type="additive", s.window=NULL, s.degree=0, t.window=NULL, t.degree=2, robust=FALSE, trend=FALSE) {
	# loess allows only an additive model
	call <- match.call()
	x <- as.ts(x)
	if (is.matrix(x)) 
        stop("only univariate series are allowed")
    # Check the type argument
	TYPES <- c("additive", "multiplicative")
		typeindex <- pmatch(type, TYPES)
		if (is.na(typeindex)) 
			stop("invalid type value")
		if (typeindex == -1) 
			stop("ambiguous type value")
		# multiplicative model not handled (yet)
		if (typeindex == 2) 
			stop("multiplicative model not handle by this function. Use deccensus() instead")
		# make sure type is fully spelled
		type <- switch(typeindex,
				"additive"="additive",
				"multiplicative"="multiplicative")
	if (type == "multiplicative")
		stop("'loess' method allows only an additive model. Use 'census' instead.")
	# create our own specs component
	specs <- list(method="loess", type=type, s.window=s.window, s.degree=s.degree, t.window=t.window, t.degree=t.degree, robust=robust, trend=trend)
	# we recuperate units from x
	units <- attr(x, "units")
	if (t.degree == 2) t.degree <- 1		# Only 0 or 1 for R
	res.stl <- stl(x, s.window=s.window, s.degree=s.degree, t.window=t.window, t.degree=t.degree, robust=robust)
	if (trend == TRUE) {
		series <- cbind(res.stl$time.series[, "trend"], res.stl$time.series[, "seasonal"], res.stl$time.series[, "remainder"])
		dimnames(series)[[2]] <- c("trend", "seasonal", "residuals")
	} else {			# residuals is trend + remainder in the additive model (otherwise, we recalculate them)
		series <- cbind(res.stl$time.series[, "trend"] + res.stl$time.series[, "remainder"], res.stl$time.series[, "seasonal"])
		dimnames(series)[[2]] <- c("deseasoned", "seasonal")
 	}
	# create our own 'tsd' structure
	res <- list(ts="series", series=series, weights=res.stl$weights, units=units, specs=specs, call=call)
	class(res) <- "tsd"		# change the class of the object to 'tsd'
	res
}
