plotVAR1data <- function(Y, lwd=1){
	######################################################################################
	# 
	# DESRIPTION:
	# Time series plot.
	#	
	# ARGUMENTS:
	# -> Y             : Three-dimensional array containing the data. The first, second and third dimensions correspond to 
	#                   covariates, time and samples, respectively. The data are assumed to centered covariate-wise.
	# -> lwd           : Numeric, specifying the width of the time course lines.
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
	if (as.character(class(lwd)) != "numeric"){ stop("Input (lwd) is of wrong class.") }
	if (length(lwd) != 1){ stop("Input (lwd) is of wrong length.") }
	if (is.na(lwd)){ stop("Input (lwd) is not a non-negative number.") }
	if (lwd < 0){ stop("Input (lwd) is not a non-negative number.") }
	
	# start plotting
	plot(y=Y[1,,1], x=1:dim(Y)[2], col="white", ylim=c(min(Y), max(Y)), xlab="time points", ylab="data")
	for (i in 1:dim(Y)[3]){
		for (j in 1:dim(Y)[1]){
			lines(1:dim(Y)[2], Y[j,,i], col=j, lty=i, lwd=lwd)
		}
	}
}
