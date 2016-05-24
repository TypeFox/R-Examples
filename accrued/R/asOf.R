asOf =  function(x, currentDate=NULL) { 
	
	## Throw an error if the argument is not of the correct class.
	if( class(x) != "accrued" )  stop("ERROR: argument is not an object of the 'accrued' class.")

	accrued_data = x		

	NROW = nrow(accrued_data[["data"]])
	NCOL = ncol(accrued_data[["data"]])
	startDate = accrued_data[["start"]]
	rowNames = as.vector(dimnames(accrued_data[["data"]])[[1]])

	if( is.null(currentDate) ) currentDate = rowNames[NROW]
	
	rowIndex = which( rowNames == currentDate )

	if( length(rowIndex) == 0 ) stop("ERROR: date is out of range.")
	else{
	
		minLaggedEncounterDate = max( c(rowIndex-NCOL+1, 1))
		finalVector = NULL		

		## If the date is more than NCOL days prior to the current date,
		## the final value is used.
		if( minLaggedEncounterDate > 1 ) finalVector = accrued_data[["final"]][1:(minLaggedEncounterDate-1)]

		laggedIndices = cbind(minLaggedEncounterDate:rowIndex, (rowIndex-minLaggedEncounterDate+1):1)
		laggedVector = accrued_data[["data"]][laggedIndices]
		
		## Returns the rowIndex-by-1 matrix of current values.
		## The column names are the dates.
		asOfMatrix = as.matrix(c(finalVector, laggedVector))
		dimnames(asOfMatrix) = list( rowNames[1:rowIndex], c("Current Values") ) 
		asOfMatrix

	}
}

	
