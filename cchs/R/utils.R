################################################################################
# Miscellaneous subroutines. 

################################################################################
# Find whether all elements of a vector are equal. 
allEqual <- function(x) {
	if (!is.vector(x)) stop("x must be a vector") 
	isTRUE(all.equal(x, rep(x[1], length(x))))
}

# Find the smallest difference between any pair of values in a vector. 
minDifference <- function(x, allowZeroResult=FALSE) {  
	differences <- c(sort(x)[-1]) - c(sort(x)[-length(x)])
	if (allowZeroResult) {
		return(min(differences))
	} else {
		return(min(differences[differences!=0]))
	}
} 

# Find whether a factor has any unused levels.
hasUnusedLevels <- function(x) {
	if (!is.factor(x)) stop("x must be a factor")
	!identical(sort(levels(x)), sort(levels(droplevels(x))))
}

# Remove NAs from a vector. 
removeNAsFromVector <- function(x) {
	if (!is.vector(x)) stop("x must be a vector")
	as.vector(na.omit(x))  # or x[!is.na(x)]
}

################################################################################
# Set the S3 class of an object. "x <- setS3class(x, className)" is the same as 
# "class(x) <- className", but it works in S as well as R. This is partly based
# on the last few lines of survival:::summary.coxph and partly based on ?is.R.
setS3class <- function(x, className) {
	if (exists("is.R") && is.function(is.R) && is.R()) {
		# this code is being run in R
		class(x) <- className 
	} else { 
		# this code is being run in S
		oldClass(x) <- className
	}  
	return(x)
}
################################################################################

