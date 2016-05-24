"tsd" <-
function (x, specs=NULL, method="loess", type=if(method == "census") "multiplicative" else "additive", lag=1, axes=1:5, order=1, times=1, sides=2, ends="fill", weights=NULL, s.window=NULL, s.degree=0, t.window=NULL, t.degree=2, robust=FALSE, trend=FALSE, xreg=NULL) {
	call <- match.call()
    x <- as.ts(x)
  	# Do we have specs?
 	if (!is.null(specs)) {
 		# Verify it is an object of the right class 'specs.tsd'
 		specs.class <- class(specs)
 		if (is.null(specs.class) || specs.class != "specs.tsd")
 			stop("specs must be a 'specs.tsd' object.")
 		# For each argument we look if it was not explicitly given (in this case, we take value in specs)
 		arg.names <- names(call)
 		if (pmatch("method", arg.names, 0) == 0)
 			method <- specs$method
 		if (pmatch("type", arg.names, 0) == 0)
 			type <- specs$type
 		if (pmatch("lag", arg.names, 0) == 0)
 			lag <- specs$lag
 		if (pmatch("axes", arg.names, 0) == 0)
 			axes <- specs$axes
 		if (pmatch("order", arg.names, 0) == 0)
 			order <- specs$order
 		if (pmatch("times", arg.names, 0) == 0)
 			times <- specs$times
 		if (pmatch("sides", arg.names, 0) == 0)
 			sides <- specs$sides
 		if (pmatch("ends", arg.names, 0) == 0)
 			ends <- specs$ends
 		if (pmatch("weights", arg.names, 0) == 0)
 			weights <- specs$weights
 		if (pmatch("s.window", arg.names, 0) == 0)
 			s.window <- specs$s.window
 		if (pmatch("s.degree", arg.names, 0) == 0)
 			s.degree <- specs$s.degree
 		if (pmatch("t.window", arg.names, 0) == 0)
 			t.window <- specs$t.window
 		if (pmatch("t.degree", arg.names, 0) == 0)
 			t.degree <- specs$t.degree
 		if (pmatch("robust", arg.names, 0) == 0)
 			robust <- specs$robust
 		if (pmatch("trend", arg.names, 0) == 0)
 			trend <- specs$trend
 		if (pmatch("xreg", arg.names, 0) == 0)
 			xreg <- specs$xreg
 	}
 	# Evaluate arguments now
 	method <- method
 	type <- type
 	lag <- lag
 	axes <- axes
 	order <- order
 	times <- times
 	sides <- sides
 	ends <- ends
 	weights <- weights
 	s.window <- s.window
 	s.degree <- s.degree
 	t.window <- t.window
 	t.degree <- t.degree
 	robust <- robust
 	trend <- trend
 	xreg <- xreg
 	nser <- ncol(x)
	if (is.null(nser)) nser <- 1		# when x is a single vector
 	# Verify the method argument
	METHODS <- c("diff", "average", "median", "evf", "reg", "census", "loess")
	methindex <- pmatch(method, METHODS)
	if (is.na(methindex)) 
		stop(paste("invalid decomposition method:", method))
	if (methindex == -1) 
		stop(paste("ambiguous decomposition method:", method))
	# make sure method is fully spelled
	method <- switch(methindex,
				"diff"="diff",
				"average"="average",
				"median"="median",
				"evf"="evf",
				"reg"="reg",
				"census"="census",
				"loess"="loess")
	# Verify the type argument
	TYPES <- c("additive", "multiplicative")
		typeindex <- pmatch(type, TYPES)
		if (is.na(typeindex)) 
			stop(paste("invalid decomposition type:", type))
		if (typeindex == -1) 
			stop(paste("ambiguous decomposition type:", type))
		# make sure type is fully spelled
		type <- switch(typeindex,
					"additive"="additive",
					"multiplicative"="multiplicative")
	# Create a specs list
	specs <- list(method=method, type=type, lag=lag, axes=axes, order=order, times=times, sides=sides, ends=ends, weights=weights, s.window=s.window, s.degree=s.degree, t.window=t.window, t.degree=t.degree, robust=robust, trend=trend, xreg=xreg)
	# Decompose each series in turn
	# If only one series, then only one call is required
	if (nser == 1) {
		# Choose the method
		res <- switch(methindex,
		    	"diff"=decdiff(x, type=type, lag=lag, order=order, ends=ends),
		    	"average"=decaverage(x, type=type, order=order, times=times, sides=sides, ends=ends, weights=weights),
				"median"=decmedian(x, type=type, order=order, times=times, ends=ends),
				"evf"=decevf(x, type=type, lag=lag, axes=axes),
				"reg"=decreg(x, xreg=xreg, type=type),
				"census"=deccensus(x, type=type, trend=trend),
				"loess"=decloess(x, type=type, s.window=s.window, s.degree=s.degree, t.window=t.window, t.degree=t.degree, robust=robust, trend=trend))
		# res is already a 'tsd' object, we just have to change specs to make sure all args are included
		res$specs <- specs
	} else {
		# If x has multiple series, we must deal with each one in turm
		res <- NULL
		# We keep the names of the series in ts
		res$ts <- dimnames(x)[[2]]
		# We initiate series
		res$series <- list(1)
		res$units <- attr(x, "units")
		res$specs <- specs
		res$call <- call
		# Calculation is performed alternatively on each column
		for (i in 1:nser) {
			res1 <- switch(methindex,
		    	"diff"=decdiff(x[,i], type=type, lag=lag, order=order, ends=ends),
				"average"=decaverage(x[,i], type=type, order=order, times=times, sides=sides, ends=ends, weights=weights),
				"median"=decmedian(x[,i], type=type, order=order, times=times, ends=ends),
				"evf"=decevf(x[,i], type=type, lag=lag, axes=axes),
				"reg"=decreg(x[,i], xreg[,i], type=type),
				"census"=deccensus(x[,i], type=type, trend=trend),
				"loess"=decloess(x[,i], type=type, s.window=s.window, s.degree=s.degree, t.window=t.window, t.degree=t.degree, robust=robust, trend=trend))
			# Add the series in res
			res$series[[i]] <- res1$series
			# If we get weights, add them too
			res$weights <- res1$weights
		}
	}
	class(res) <- "tsd"
	res
}
