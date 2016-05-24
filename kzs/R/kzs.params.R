kzs.params <- function(x, dimension) {
	smooth <- numeric(dimension)
	desc <- character(dimension)
	scale <- numeric(dimension)
	ddesc <- character(dimension)
	for (i in 1:dimension) {
		smooth[i] <- max(x[,i]) - min(x[,i])
		sx <- sort(x[,i])
		dx <- diff(sx)
		scale[i] <- min(dx[dx > 0])
		desc[i] <- paste("For x",i,", 'smooth' must be a positive real number much less than ", smooth[i], sep = "")
		ddesc[i] <- paste("For x",i,", 'scale' must be a positive real number less than ", scale[i], sep = "")
	}	
	lst <- list(smooth = desc, scale = ddesc)
	return(lst)
}