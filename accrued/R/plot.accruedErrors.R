##########################################################################################
##########################################################################################

plot.accruedErrors = function(  x, 
			    				withSmoothing = FALSE,
			    				quantiles = c(0.1, 0.5, 0.9), 
			    				quantileColors = switch(1+is.null(quantiles), rainbow(length(quantiles)), NULL), ...)
{

	if( class(x) != "accruedErrors")  stop("ERROR: the first argument is not of 'accruedErrors' class.")
	
	STACKED = x
	LAG_QUANTILES = errorQuantileSummary(x, quantiles=quantiles)

	################################################################################################################
	## Plot the quantiles.
	################################################################################################################

	Y_MAX = max( STACKED[,"Error"], na.rm=T ) 
	Y_MIN = min( STACKED[,"Error"], na.rm=T ) 
	if( Y_MAX < 0 ) Y_MAX = 0.5
	if( Y_MIN > 0 ) Y_MIN = -0.5
	X_MIN = 0
	X_MAX = max( STACKED[,"Lag"], na.rm=T ) 
	plot( c(X_MIN,X_MAX), c(Y_MIN,Y_MAX), 
		  xlim=c(X_MIN,X_MAX), ylim=c(Y_MIN,Y_MAX), 
		  xlab="Lag", ylab="Error",
		  type = 'n',  pch=16, cex=0.3, 
		  col="dimgrey",
		  axes=F,  
		  main="", ... )
	abline( a=0,b=0, lwd=0.5, col="gray" )
	axis(1)
	axis(2, las=2)
	points( jitter(STACKED[,"Lag"]), STACKED[,"Error"] , pch=16, cex= 0.3, col="dimgrey")
	if (!is.null(quantiles))  
		for( q in 1:length(quantiles) )  
			if( withSmoothing ) {
				if( sum(is.na(as.vector(LAG_QUANTILES))) == 0 ) 
					lines( smooth.spline(as.numeric(colnames(LAG_QUANTILES)), LAG_QUANTILES[q,] ), col=quantileColors[q] ) 
			} else lines( as.numeric(colnames(LAG_QUANTILES)), LAG_QUANTILES[q,], col=quantileColors[q] )

	invisible(LAG_QUANTILES)
	
} # END quantile function 
