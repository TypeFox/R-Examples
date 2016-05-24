###########################
## Transform dist object to a matrix using fast and efficiant C code
###########################

dist2matrix <- function(dist) {
	if (inherits(dist, "dist")) {
		return(.Call(TMR_dist2matrix, as.double(dist), attr(dist, "Size")))
	}
	else if (is.matrix(dist)) {
		return(dist)
	}
	stop("dist should be a matrix or a \"dist\" object")
}