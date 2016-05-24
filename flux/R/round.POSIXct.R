round.POSIXct <- function(x, ...){
	x <- as.POSIXlt(x)
	return(round.POSIXlt(x, ...))
}