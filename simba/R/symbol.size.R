"symbol.size" <- function(x, cex.min = 0.2, cex.max = 5){
	r.y <- range(x, na.rm=TRUE)
	size <- cex.min + ((x - r.y[1]) * (cex.max - cex.min))/(r.y[2] - r.y[1])
	size[is.na(size)] <- cex.min
	size
}