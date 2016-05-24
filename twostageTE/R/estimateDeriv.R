estimateDeriv <-
function (explanatory, response, d_0, sigmaSq) {

	deriv_estimateHelper <-
	function(explanatory, response, d_0, sigmaSq) {
		n <- length(response)
		p <- 5 # 5th order needed 	
		# X is the design mtx
		X <- matrix(0, n, p)
	#	X[,1] <- 1
		for (i in 1:p) {
			X[,i] <- (explanatory - d_0)^i
		}
	#	Y <- response
	#	# Now construct beta_hat
	#	m1 <- crossprod(X,X)
	#	m2 <- crossprod(X,Y)
	#	# Compute the inverse
	#	m1_inv <- solve(m1)
	#	# Compute the coefficient vector of the weighted LSE problem.
	#	beta_hat <- m1_inv %*% m2
		beta_hat <- lm(response ~ 0+X)$coef
		h <- 0
		for (i in (p-1):(p+1)) {
			j <- i - p + 2 
			h <- h + beta_hat[i-1]*factorial(j)*d_0^(j-1)
		}
#		sigmaSq <- estimateSigmaSq(explanatory, response)$sigmaSq
		# Return the optimal bandwidth.
		return(2.275 * (sigmaSq / h^ 2.0) ^ (1.0/7.0) * n ^ (-1.0/7.0))
	}


	n <- length(response)
	p <- 2 # need only quadratic term since we're interested in the first derivative
	X <- matrix(0, n, p) 	# X is the design mtx
	X[,1] <- (explanatory - d_0)
	X[,2] <- (explanatory - d_0)^2

	bw_opt <- deriv_estimateHelper(explanatory, response, d_0, sigmaSq)
	# Construct the weight matrix with Epanechnikov kernel
#	W <- matrix(0,n,n)
#	for (i in 1:n)   #### SHOULD VECTORIZE THIS TO RID FOR LOOP
#		W[i,i] <- 0.75 * max(1.0-( (explanatory[i] - d_0)/ bw_opt )^2, 0) /bw_opt 
#print(W)
#	Y <- response
#	# Now construct beta_hat
#	tmp <- crossprod(X,W)
#	m1 <- tmp %*% X
#	m2 <- tmp %*% Y
#	# Compute the inverse
#	m1_inv <- solve(m1)
#	# Compute the coefficient vector of the weighted LSE problem.
#	beta_hat <- m1_inv %*% m2
#print(beta_hat)
	W <- 0.75 /bw_opt  * sapply(1.0-( (explanatory - d_0)/ bw_opt )^2, max,0)
	while (sum(W>1) <= 1 & bw_opt <= max(explanatory) - min(explanatory)) {
		bw_opt <- bw_opt*2
		W <- 0.75 /bw_opt  * sapply(1.0-( (explanatory - d_0)/ bw_opt )^2, max,0)	
	}
	beta_hat <- lm(response ~ 0+X, weights=W)$coef
 	while (beta_hat[1] <= 0 & bw_opt <= max(explanatory) - min(explanatory)) {
		bw_opt <- bw_opt*2
		W <- 0.75 /bw_opt  * sapply(1.0-( (explanatory - d_0)/ bw_opt )^2, max,0)	
		beta_hat <- lm(response ~ 0+X, weights=W)$coef
	}
 	if (beta_hat[1] <= 0) {warning("Negative derivative has been estimated. Methods within twostageTE assume underlying non-decreasing functions. Consider transforming your data so that it is non-decreasing." , call.=FALSE)}
	return (beta_hat[1]) 
}
