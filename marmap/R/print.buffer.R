print.buffer <- function(x, ...){
	
	if (!is(x,"buffer")) stop("x must be an object of class buffer")
	
	cat("Object of class 'buffer' containing:\n")
	cat("1. a 'bathy' object filled with NAs except inside the buffer\n")
	cat("2. a two column data.frame locating the center of the buffer\n")
	cat("3. the radius of the buffer in the same unit as the 'bathy' object\n")
	cat("4. the radius of the buffer in kilometers\n")
	cat("\n")
	cat("Summary of the 'bathy' object\n")
	print(summary(x[[1]]))
	cat("\n")
	cat("Coordinates of the center of the buffer\n")
	print(x[[2]])
	cat("\n")
	cat("Radius of the buffer\n")
	print(x[[3]])
	cat("\n")
	cat("Radius (in kilometers) of the buffer\n")
	print(x[[4]])
}