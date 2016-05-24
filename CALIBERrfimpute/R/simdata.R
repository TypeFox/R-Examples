simdata <-
function(n=2000, mymean=rep(0,4), mysigma=matrix(
	c(1  , 0.2, 0.1,-0.7,
	  0.2, 1  , 0.3, 0.1,
	  0.1, 0.3, 1  , 0.2,
	 -0.7, 0.1, 0.2, 1), byrow=TRUE, nrow=4, ncol=4),
	residsd=1, x2binary=FALSE){
	# Returns a simulated dataset. Predictors are drawn from a multivariate
	# normal distribution with mean and covariance sigma, with residual
	# variance residsd

	# Covariance matrix: 1    0.2  0.1 -0.7
	#                    0.2  1    0.3  0.1
	#                    0.1  0.3  1    0.2
	#                   -0.7  0.1  0.2  1
	
	# output is calculated based on the predictors
	# two of interest, two auxiliary
	out <- rmvnorm(n, mymean, mysigma)
	if (x2binary==TRUE){
		# Convert x2 to a random draw from the logistic of x2
		out[,2] <- rbinom(n, 1, exp(out[,2])/(1+exp(out[,2])))		
	}
	# y is the sum of the first 3 x variables (i.e. true coefficients are 1)
	# add a random error to the output
	out <- cbind(rnorm(n, out[,1] + out[,2] + out[,3], residsd), out)
	dimnames(out)[[2]] <- c('y', paste('x', 1:4, sep=''))
	out <- data.frame(out)
	if (x2binary==TRUE){
		out$x2 <- as.factor(out$x2 + 1)
	}
	return(out)
}
