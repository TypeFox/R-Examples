data.accrued <- function( data,   
						 start = NULL, 
						 final = NULL ) {

	####################################################################################################
	## If the dimensions are too small, exit.
	if( dim(data)[[2]] == 0 ) stop("ERROR: in 'accrued' data matrix is empty.")

	###############################
	## INITIALIZING COLUMN NAMES ##
	###############################
	lag.r = ncol(data) - 1
	colnames(data) = paste( c( 0:lag.r ) )

	###############################
	## INITIALIZING FINAL COUNTS ##
	####################################################################################################	
	## If "final" the vector of final numbers (counts) isn't in the input, then 
	## the final values to be the right-most final column of data.
	####################################################################################################	
	final.r = final
	if( is.null(final) ) {
		FINAL_LABEL = paste( lag.r, sep="" )	 
		final.r = data[, FINAL_LABEL] 
	}
	
	##########################################
	## INITIALIZING "start" ENCOUNTER DATES ##
	##########################################
	start.r = start
	end.r = NULL
	if( is.null(start) ) {
		start.r = 1
		end.r = nrow(data)
	} else {
		start.r = as.Date(start)
		end.r = as.Date( start.r + nrow(data) - 1 )
	}

	times.r = seq( start.r, end.r, by=1 )
	rownames(data) = as.character( times.r )			
			
	# start.r is either a date or the integer 1.
	result = list( final=final.r, data=data, start=start.r )	
	class(result) = 'accrued'
	result
}

####################
## Print function ##
####################
print.accrued <- function(x, ...)  {
	## Throw an error if the argument is not of the correct class.
	if( class(x) != "accrued" )  stop("ERROR: argument is not an object of the 'data.accrued' class.")
	out = cbind( x$data, x$final )
	colnames(out) = c( colnames(x$data), "final" )
	rownames(out) = rownames(x$data)
	print(out, ...)
}





