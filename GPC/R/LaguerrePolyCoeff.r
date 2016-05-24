LaguerrePolyCoeff <- function(degree,alpha){
	if (degree == 0){
		value <- c(1.0, 0.0, 0.0)
	} else {
		a <- -1.0/degree
		b <- (2.0*degree-1.0+alpha)/degree
		c <- -(degree-1.0+alpha)/degree
		value <- c(a,b,c)
	}
	return (value)
}
