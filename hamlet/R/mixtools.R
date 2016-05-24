#
# Tools for the analysis of mixed data
#

# Apply an R-function to the numeric columns of a data.frame
# Useful for such functions as:
# scale - scaling numeric columns (x - mean / sd)
# cov - computing covariance matrix
# cor - computing correlation matrix
mix.fun <- function(
	# Data frame or matrix to scale
	x,
	# Function to be used
	FUN = scale,
	# Additional parameters passed to the R function FUN
	...
){
	# Checking which columns are numeric
	isnum <- lapply(x, is.numeric)
	# Performing scaling for columns that are numeric
	nums <- which(isnum==TRUE)
	# Issue a warning if nothing is done
	if(length(nums)==0){
		warning("No numeric columns detected for computing covariance matrix!")
	}
	return <- FUN(x[,nums], ...)
	# Returning the matrix with numeric columns changed
	return
}


# A function that codes categorical variables in a dataset into binary variables. This is done in the following manner:
# x = {red, green, blue, green} --> x_new = {{1,0,0}, {0,1,0}, {0,0,1}, {0,1,0}} where the dimensions in x_new are is_red, is_green and is_blue
mix.binary <- function(
	# Data.frame
	x
){
	# Checking which columns are numeric
	isnum <- lapply(x, is.numeric)
	# Performing coding for columns that are non-numeric
	nonnums <- which(isnum==FALSE)
	# Issue a warning if nothing is done
	if(length(nonnums)==0){
		warning("Only numeric columns detected, unable to process binary coding of categorical variables!")
	}
	xnew <- x
	# Looping through categorical variables in the original data matrix
	for(i in 1:length(nonnums)){
		uniqs <- unique(x[,nonnums[i]])
		newcols <- paste(colnames(x)[nonnums[i]],"_", uniqs, sep="")
		#print(newcols)
		bins <- matrix(0, nrow=nrow(x), ncol=length(uniqs))
		colnames(bins) <- newcols
		# Looping through observations
		for(j in 1:nrow(x)){
			ind <- which(x[j,nonnums[i]]==uniqs)
			bins[j,ind] <- 1
		}
		xnew <- cbind(xnew, bins)	
	}
	# Removing the old, non-binary coded fields
	xnew <- xnew[,-c(nonnums)]
	
	# Returning the matrix with the newly binary coded fields
	xnew
}

