plot.buffer <- function(x, outline = TRUE, add = TRUE, ...) {

	buffer <- x  # S3 compatibility

	if (!is(buffer,'buffer')) stop("'buffer' must be an object of class bathy as produced by create.buffer()")

	ll <- list(...)

	if (outline) {
		symbols(buffer[[2]], circles = buffer[[3]], add = add, inches = F, ...)
	} else {
		plot.bathy(buffer[[1]], add=add, ...)
	}

}