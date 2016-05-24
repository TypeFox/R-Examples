centerVAR1data <- function(Y){
	######################################################################################
	#
	# DESCRIPTION:
	# within-individual, covariate-wise centering of the data.
	#
	# ARGUMENTS:
	# -> Y             : Three-dimensional array containing the data. The first, second and third dimensions correspond to 
	#                   covariates, time and samples, respectively. The data are assumed to centered covariate-wise.
	#
	# DEPENDENCIES:
	# ...
	#
	# NOTES:
	# ...
	#
	######################################################################################

	# input checks
	if (as.character(class(Y)) != "array"){ stop("Input (Y) is of wrong class.") }
	if (length(dim(Y)) != 3){ stop("Input (Y) is of wrong dimensions: either covariate, time or sample dimension is missing.") }

	# centering
	for (i in 1:dim(Y)[3]){
	    Y[1:dim(Y)[1], 1:dim(Y)[2],i] <- sweep(Y[,,i,drop=TRUE], 1, apply(Y[,,i,drop=TRUE], 1, mean, na.rm=TRUE))
	}
	return(Y)
}


