`sinc` <-
function(X)
{
	# copy X to Y, so Y has same format as X
	Y <- X
	
	# find zero elements in X
	Z <- X == 0
	Y[Z] <- 1 # sinc(0) = 1
	Y[!Z] <- sin(pi*X[!Z]) / (pi*X[!Z]) # sinc(x) = sin(x)/x for x != 0
	
	return(Y)
}

