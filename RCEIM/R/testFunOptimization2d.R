testFunOptimization2d <- function(x) {
# The function to be optimized
# boundaries <- matrix(nrow=2, ncol=2)
# boundaries[1,] <- c(-10, 10)
# boundaries[2,] <- c(-10, 10)
# 
	return( ((x[1]-4)^2 + (x[2]+2)^2)/50 - ((x[1]+2)^2 + (x[2]+4)^2)/90 
		-exp(-((x[1]-2)^2)) - 0.9*exp(-((x[2]+2)^2)) - 0.5*sin(x[1]*8) - 0.25*cos(x[2]*2) 
		+0.25*sin(x[1]*x[2]/2) + 0.5*cos(x[1]/2.5*x[2]))
}