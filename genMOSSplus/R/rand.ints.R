rand.ints <- function(n, min=1, max=100) {
#
# Randomly generates EXACTLY n integers in the range from min to max, inclusive.
#

total <- max - min + 1

if(n <= 0) {
	print(paste("Warning: the number of unique integers required is too small: ", n, sep=""))
	return(numeric(0))
}
if(n > (total)) {
	print(paste("Warning: the number of unique integers required: ", n, " exceeds the size of the interval: [", min, ", ", max, "].", sep=""))
	return(seq(from=min, to=max, by=1))
}


# It will be faster to sample fewer random numbers.
# Thus if more than 50% of the numbers are required within the [min,max] range,
# we will find the opposite set (the one to throw out).
# This should reduce the number of times the while loop will get executed.

p <- n/total*100  

if(p > 50) {
	
	vals <- helper.rand.ints(total - n, min, max)
	vals <- setdiff(min:max, vals)
} else {
	vals <- helper.rand.ints(n, min, max)
}

return(vals)
}


# Helper function that actually does the sampling of n random values.
# Would be best if the number of values to generate, n, is smaller
# than half of the [min, max] range. 
# Also n must be greater than zero, and less then total number of integers within the range.
# Returns the array of length n containing random integers in [min,max] range.

helper.rand.ints <- function(n, min, max) {

vals <- unique(round(runif(n, min=min, max=max)))

while(length(vals) < n) {
	new.val <- round(runif(1, min=min, max=max))
	val.dne <- is.na(match(new.val, vals))
	if(val.dne==TRUE)
		vals <- union(vals, new.val)
}

return(vals)

}



