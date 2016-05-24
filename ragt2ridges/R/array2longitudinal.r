array2longitudinal <- function(Y, keepMissings=TRUE){
	######################################################################################
	#
	# DESCRIPTION:
	# Converts a 3-dim array (containing time-series data of multiple individuals) to an
	# object of which can directly by converted to the 'longitudinal' class.
	#
	# ARGUMENTS:
	# -> Y			: Three-dimensional array containing the data. The first, second 
	# 			  and third dimensions correspond to covariates, time and 
	#			  samples, respectively.
	# -> keepMissings	: The 'array'-format assumes a balanced lay out of the 
	# 			  time-course experiment. The experiment may have failed for 
	# 			  some design and no data is available. In the 'longitudinal'-
	# 			  format these design points may be left out. This 'logical'
	# 		          indicates whether they should be kept in (or not).
	#
	#
	# NOTES:
	# ...
	#
	######################################################################################

	rNames <- paste(sort(rep(1:dim(Y)[2], dim(Y)[3])), rep(1:dim(Y)[3], dim(Y)[2]), sep="-")
	Ymatrix <- matrix(NA, ncol=dim(Y)[1], nrow=dim(Y)[2] * dim(Y)[3])
	for (t in 1:dim(Y)[2]){
		Ymatrix[(t-1)*dim(Y)[3] + 1:dim(Y)[3], ] <- t(Y[,t,])
	}
	rownames(Ymatrix) <- rNames
	attr(Ymatrix, "time") <- 1:dim(Y)[2]
    	attr(Ymatrix, "repeats") <- rep(dim(Y)[3], dim(Y)[2])
	return(Ymatrix)
}

