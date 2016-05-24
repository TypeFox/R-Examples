
###############################################################################################################
## This function creates stacked upload data. Its argument is
## an object of type 'data.accrued' and it returns a matrix in which "lag" is treated
## as a variable. It is used for the sparkline plots ('plot.accrued'), the barcode plots, and upload plots.
###############################################################################################################


stackedUploadData = function(x) {

	## Throw an error if the argument is not of the correct class.
	if( class(x) != "accrued" )  stop("ERROR: argument is not an object of the 'accrued' class.")

	accrued_data = x

	data = accrued_data[["data"]]
	MAX_LAGS = ncol(data) - 1
	NCOL = ncol(data)
	NROW = nrow(data)

	## Data structure to be returned.
	STACKED = matrix(  NA, nrow=NROW*NCOL, ncol=5, 
					   dimnames=list( 1:(NROW*NCOL), 	
					   c( "EncounterDate", "UploadDate", "Lag", "Counts", "NumberAdded") )  )

	## Populate encounter dates.
	STACKED[,"EncounterDate"] = rep( 1:NROW, times=NCOL )

	## Populate lags.
	STACKED[,"Lag"] = rep(0:MAX_LAGS, rep(NROW, MAX_LAGS + 1))
	
	## Populate upload dates.
	STACKED[,"UploadDate"] = STACKED[,"EncounterDate"] + STACKED[,"Lag"] 

	## Populate cumulative counts for that (encounter date, lag) combination.
	## Generate all possible combinations of (encounter date, lag).
	## Save these indices.
	index = cbind( STACKED[,"EncounterDate"], STACKED[,"Lag"] + 1 )

	# Populate the counts column using the original data.
	STACKED[,"Counts"] = data[index]
	
	## Indicator MATRIX of missing data (1 means data is present, 0 means data is missing).
	isNAdata = 0 + !is.na(data)
	
	## Find the first column index for which there IS data.
	## This returns a column VECTOR or length NROW that contains the index of the first column in that row
	## for which there IS a data entry.
	col.first.not.na = apply(isNAdata, 1, which.max)
	
	## For each row, find the last column index for which there is not data.
	## This returns a vector of length NROW.
	col.last.na = col.first.not.na - 1
	
	## For every column for which no data are missing, set the last column index for which
	## there is not data to be NA (this would be any index at 0).  
	## This returns a vector of length NROW.
	## The only change from the previous is that all zeroes become NAs.
	col.last.na[col.last.na < 1] = NA

	## Identify the column indices of rows for which data ARE missing (meaning, the col.last.na has a nontrivial index).
	## This creates a vector of indices including only every ROW index for which col.last.na is NOT NA.
	## Therefore it is a vector of length less than or equal to NROW. The only way for it to have
	## length NROW is if no entries of col.last.na are NA. This is equivalent to every row missing at least some data.
	## And, if that's the case, then adjust.rows will just be the vector 1:NROW.
	adjust.rows = which( !is.na(col.last.na) )
	
	## Create a temporary copy of the data.
	temp_data = data
	
	## Only for the indices stored in 'adjust.rows'
	## Set the the right-most occuring "NA" in that location equal to "0".
	## The preceding NAs need not be reset since they won't be used to recover the number added. 
	## On that upload date.
	for( i in adjust.rows )  temp_data[i, col.last.na[i]] = 0
	

	## This is the line that fails if there is only one column to begin with.		
	## The "else" of the if-else returns the number of counts added.
	## The first column has nothing added from the previous since it's the first upload date for that encounter date.
	## Then, to determine the number of counts added, the temporary data with columns of lag 1 to the end has subtracted
	## from it the temporary data with columns of lag 0 to end-1 (or MAX_LAGS). 
	## Matrix subtraction is faster entry-wise subtraction using for loops.
	data.lag = NULL
	if( MAX_LAGS == 0 ) numberAddedMatrix = data[,1, drop=FALSE]
	else numberAddedMatrix = cbind(  data[,1],  temp_data[, 2:NCOL] - temp_data[,1:MAX_LAGS]  )

	## This assigns the column "NumberAdded" the stacked vector "numberAddedMatrix", which stackes left-column to right-column.
	STACKED[,"NumberAdded"] = numberAddedMatrix[index]

	as.data.frame(STACKED)

}
