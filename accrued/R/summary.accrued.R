
## This creates the accrued summary.
summary.accrued = function(object, ...) {

	## Throw an error if the argument is not of the correct class.
	if( class(object) != "accrued" )  stop("ERROR: argument is not an object of the 'accrued' class.")
		
	x = object
	lags = ncol(x$data)

	# An anonymous function is created in the 3rd argument of the "apply" function call below.
	# "y" is a generic object (whatever "apply" passes to the function) and returns the proportion
	# of non-missing data. (!is.na(y) returns T (which is then coerced to 1) if data are present and F otherwise)

	## "upload.prop" is created as follows. For each lag, an encounter date is assigned a 1 if data for that encounter date,
	## are uploaded, and a 0 otherwise. 
	## The percentage of such encounter dates is returned for each lag.
	upload.prop = apply( cbind(x$data, x$final), 2, function(y) mean(!is.na(y)) )


	## "mean.prop" is created as follows. For each lag, the proportion of counts for reach upload date, 
	## relative to the final counts, is computed. Then, the mean over encounter dates is taken, 
	## for each lag, and returned.
	mean.prop   = apply( cbind(x$data, x$final), 2, function(y) mean(y/x$final, na.rm=TRUE) )


	## "q1.prop" is created as follows. For each lag, the proportion of counts for reach upload date, 
	## relative to the final counts, is computed. Then, the 1st quartile over encounter dates is taken, 
	## for each lag, and returned.
	q1.prop   = apply( cbind(x$data, x$final), 2, function(y) quantile(y/x$final, probs=0.25,na.rm=TRUE)  )

	## "q2.prop" is created as follows. For each lag, the proportion of counts for reach upload date, 
	## relative to the final counts, is computed. Then, the 2nd quartile over encounter dates is taken, 
	## for each lag, and returned.
	q2.prop   = apply( cbind(x$data, x$final), 2, function(y) quantile(y/x$final, probs=0.5,na.rm=TRUE)  )

	## "q3.prop" is created as follows. For each lag, the proportion of counts for reach upload date, 
	## relative to the final counts, is computed. Then, the 3rd quartile over encounter dates is taken, 
	## for each lag, and returned.
	q3.prop   = apply( cbind(x$data, x$final), 2, function(y) quantile(y/x$final, probs=0.75,na.rm=TRUE)  )

	## "mean.tot" is created as follows. The mean number of counts for any encounter date is computed, for
	## each lag.
	mean.tot = apply( cbind(x$data, x$final), 2, mean, na.rm = TRUE )
	
	result = data.frame( 	upload.prop=upload.prop, 
							mean.prop=mean.prop, 
							mean.total=mean.tot, 
							q1.prop=q1.prop, 
							q2.prop=q2.prop, 
							q3.prop=q3.prop    )
	row.names(result) = c(0:(lags-1), "final")
	class(result) = 'summary.accrued'
	result
}



## This prints the accrued summary.
print.summary.accrued = function(x, ...) {

	if( class(x) != "summary.accrued" )  stop("ERROR: argument is not an object of the 'summary.accrued' class.")

	upload.prop = x$upload.prop
	mean.prop = x$mean.prop
	mean.tot = x$mean.tot
	q1.prop = x$q1.prop
	q2.prop = x$q2.prop
	q3.prop = x$q3.prop
	
	lags = length(upload.prop)-1
	out = cbind( 	format(round(upload.prop, 2)), 
					format(round(mean.prop, 2)), 
					format(round(mean.tot, 1)),
					format(round(q1.prop, 2)), 
					format(round(q2.prop, 2)), 
					format(round(q3.prop, 2)) )
	out = rbind( c("","","", "","",""), out )
	row.names(out) = c(" Lag", paste(rep('', lags), as.character(0:(lags - 1))), " final")
	dimnames(out)[[2]] = c( "Upload Proportion", "Mean Proportion", "Mean Count", "Quartile 1", "Median", "Quartile 3")
	
	cat( paste("\nSummary of accrued data object with", nrow(x$data), "time points.\n\n") )
	print.default( out, print.gap=2, quote=FALSE )
}





# This plots the summary.
plot.summary.accrued = function(x, ...) {
	if( class(x) != "summary.accrued" )  stop("ERROR: argument is not an object of the 'summary.accrued' class.")
	lags = length(x$mean.prop)
	askValue = par()$ask
	par(ask=TRUE)
	plot( 0:(lags-2), x$mean.prop[-lags], ylim=c(0,1), type='l', xlab='lag', ylab='proportion of counts received', 
		main="mean proportion vs lag", axes=F,...)
	axis(1)
	axis(2,las=2)
	plot( 0:(lags-2), x$q2.prop[-lags], ylim=c(0,1), type='l', xlab='lag', ylab='proportion of counts received',
		main="1st, 2nd and 3rd quartiles vs lag", axes=F,...)
	axis(1)
	axis(2,las=2)
	lines(0:(lags-2), x$q1.prop[-lags], col="forestgreen")
	lines(0:(lags-2), x$q3.prop[-lags], col="dodgerblue")
	par(ask = askValue)
}




