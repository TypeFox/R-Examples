result.extract.sub <- function(mask.grid, values, gk4.x, gk4.y, 
				outliers, silent=FALSE, withOutliers=FALSE){
	#prepare output-message
	if (!silent){ msg <- "" }

	# prepare vectors
	extracted.values <- rep(NA, times=length(mask.grid$alt))

	# iterate over all values
	for (i in 1:length(values)){
		# get gridcellnumber out of gk-coordinates
		gridcellnumber <- data.coordinates2gridcellnumber(grid=mask.grid, 
					x=gk4.x[i], y=gk4.y[i])
		
		# substitute values in prepared vectors
		if ((withOutliers)||(outliers[i]==0)){
			extracted.values[gridcellnumber] <- values[i]
		}
		# output-message
		if (!silent){
			cat(rep("\b", nchar(msg)),sep="")
			msg <- paste(i," of ", 
				length(values), " Datasets done!",sep="")
			cat(msg,sep="")
		}
	}
	if (!silent){ cat("\n") }
		
	# create data.frame
	result.values <- data.frame(values=extracted.values, 
			x=mask.grid$x, y=mask.grid$y)

	return(result.values)
}
