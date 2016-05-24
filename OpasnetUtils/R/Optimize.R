####################
# Minimize
##########################
# Inputs:
#   data - data.frame or ovariable
#   indices - character vector giving indices to be retained, all others will be minimized over
# Output:
#   data.frame which is a slice from the original
######################

Minimize <- function(data, indices){
	if (class(data) == "data.frame") data <- Ovariable(output = data)
	wm <- tapply(result(data), data@output[indices], which.min)
	
	finder <- function(vec, wm) {
		vec[wm]
	}
	
	row <- tapply(1:nrow(data@output), data@output[indices], list)
	
	out <- NULL
	for (i in 1:length(row)) {
		out <- c(out, row[[i]][wm[i]])
	}
	
	return(data@output[out,])
}

Optimize <- function(...){
	return(Minimize(...))
}