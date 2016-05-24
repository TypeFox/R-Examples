#' Internal method for plotting.
#' Position along the left side segement
#' 
#' @param d the distance
#' @export
segment1 <- function(d) {
	x <- cos((60 * pi) / 180) * d
	y <- x * (sqrt(.75) / .5)
	return(c(x=x, y=y))
}

#' Internal method for plotting.
#' Position along the right side segement
#' 
#' @param d the distance
#' @export
segment2 <- function(d) {
	x <- 1 - (cos((60 * pi) / 180) * d)
	y <- -(sqrt(.75) / .5) * x + sqrt(.75) / .5
	return(c(x=x, y=y))
}

#' Internal method for plotting.
#' Finds a point d distance from x, y
#' 
#' @param x x coordinate
#' @param y y coordinate
#' @param d the distance
#' @export
perpPt <- function(x, y, d=.05) {
	x2 <- cos( (30*pi) / 180) * d
	y2 <- sin( (30*pi) / 180) * d
	return(c(x=x2, y=y2))
}
