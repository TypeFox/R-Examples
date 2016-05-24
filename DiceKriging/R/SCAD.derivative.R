SCAD.derivative <- function(x, lambda) {
	
	x <- abs(x)
	a <- 3.7
	
	u <- (x<=lambda)
	penalty.derivative <- lambda * u + (pmax(a*lambda-x,0)/(a-1)) * (1-u)
	return(penalty.derivative)
}