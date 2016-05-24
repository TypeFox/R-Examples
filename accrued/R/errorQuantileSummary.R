
errorQuantileSummary = function( x, quantiles = c(0.1, 0.5, 0.9), ... ) {

	if( class(x) != "accruedErrors")  stop("ERROR: the first argument is not of 'accruedErrors' class.")

	STACKED = x

	LAGGED_QUANTILES = matrix(NA, nrow=length(quantiles),  ncol=(1 + max(STACKED[,"Lag"])), 
				           		  dimnames=list( quantiles, 0:max(STACKED[,"Lag"]) )  )
				       
	for ( INDEX in 1:length(quantiles) ) {
		AGGREGATED_QUANTILES =  aggregate(  STACKED[,"Error"], 
		  									list(STACKED[,"Lag"]), 
		  									quantile, 
		  									quantiles[INDEX], 
		  									na.rm = TRUE )
		LAGGED_QUANTILES[INDEX, as.character(AGGREGATED_QUANTILES[,1])] = AGGREGATED_QUANTILES[,2]
	}
	
	LAGGED_QUANTILES

}


