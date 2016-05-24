testFunOptimization <- function(x) {
# The function to be optimized
	return( exp(-((x-2)^2)) + 0.9*exp(-((x+2)^2)) + 0.5*sin(x*8) + 0.25*cos(x*2) )
}