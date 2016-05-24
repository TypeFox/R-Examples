SCAD <- function(x, lambda) {
	
	x <- abs(x)
	a <- 3.7
	
	u <- (x<=lambda)
	penalty <- lambda*x*u + (-0.5*pmax(a*lambda-x,0)^2/(a-1)  + 0.5*(a+1)*lambda^2)*(1-u)
	return(penalty)
}