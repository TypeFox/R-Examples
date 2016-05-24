"cluster" <- function(x, dist="Gamma")
 {
	if (!(dist %in% c("Gamma","LogNormal"))) { stop("Only 'Gamma' and 'LogNormal' distributions for frailties are allowed") }
	attr(x,"type") <- dist
	x
 }