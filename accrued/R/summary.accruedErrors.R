
summary.accruedErrors = function(object, ...) {

	if( class(object) != "accruedErrors")  stop("ERROR: the first argument is not of 'accruedErrors' class.")

	LAGGED_QUANTILES = errorQuantileSummary(object)
	
	LAGGED_QUANTILES

}


