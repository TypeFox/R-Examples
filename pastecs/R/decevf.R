"decevf" <-
function(x, type="additive", lag=5, axes=1:2) {
	x <- as.ts(x)
	if (is.matrix(x) && ncol(x) != 1) 
		stop("only univariate series are allowed")
	if (!is.numeric(axes) || any(axes <= 0))
		stop("axes must be a vector of positive numbers (ex 1:3)")
	if (!is.numeric(lag) || lag <= 0 || lag < max(axes))
		stop("lag must be a positive number higher or equal to axes max value")
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
	specs <- list(method="evf", type=type, lag=lag, axes=axes)
	# we recuperate units from x
	units <- attr(x, "units")
	# perform filtering
	# Create the matrix with lagged series from 0 to lag
	xlagmat <- embed(x, lag)
	# Perform a pca decomposition of this matrix
	x.pca <- princomp(xlagmat)
	# Rotated vectors are obtained by:
	# sweep(x, 2, x.pca$center) %*% x.pca$loadings == predict(x.pca)
	# original vectors are recalculated with:
	# sweep(predict(x.pca) %*% solve(x.pca$loadings, 2, x.pca$center, FUN="+")
	# for evf, we just keep some of the components in solve(x.pca$loadings)
	invloadings <- solve(x.pca$loadings)		# inverse of loadings matrix, i.e., eigenvectors
	settonul <- is.na(match(1:lag, axes))
	invloadings[settonul,] <- 0					# those are the component we drop
	xlagmat.recalc <- sweep(predict(x.pca) %*% invloadings, 2, x.pca$center, FUN="+")
	# Then we need to take the mean for diagonals to calculated filtered values of initial series
	xmat.recalc <- matrix(NA, nrow= length(x), ncol=lag)
	n <- nrow(xlagmat.recalc)
	for (i in 1:lag)
		xmat.recalc[1:n+(lag-i), i] <- xlagmat.recalc[,i]
	# perform column means to get filtered time series
	filtered <- apply(xmat.recalc, 1, mean, na.rm=TRUE)
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
