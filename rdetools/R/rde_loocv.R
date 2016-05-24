`rde_loocv` <-
function(K, y, est_y = FALSE, alldim = FALSE, est_noise = FALSE, regression = FALSE, nmse = TRUE, dim_rest = 0.5)
{
	# test if K is a matrix
	if(!is.matrix(K))
	{
		print("K must be a kernel matrix")
		return()
	}
	# test if K is a square matrix
	if(nrow(K) != ncol(K))
	{
		print("Kernel matrix K must be square matrix")
		return()
	}
	# test if y has same dimension as K
	if(length(y) != nrow(K))
	{
		print("Length of label vector y must be equal to the number of rows (and columns) of K")
		return()
	}
	
	n <- length(y) # number of samples
	
	# test if y is a column vector
	if(sum(nrow(y)) != length(y))
	{
		y <- matrix(y, n, 1) # y as a column vector
	}

	# calculate eigendecomposition of K
	E <- eigen(K, symmetric = TRUE)
	U <- E$vectors

	# start values for iteration
	s <- matrix(rep(0, n), n)
	z <- t(U) %*% y # kernel pca coefficients
	errors <- rep(0, n)

	Yh <- denoiselabels(0, U, z)
	# iterate
	for(i in 1:n)
	{
		u <- U[, i, drop = FALSE]
		s <- s + u*u
		s[(Yh[,i,drop = FALSE] - y) == 0] <- 0 # avoid zero by zero division
		errors[i] <- sum(((Yh[,i,drop = FALSE] - y) / (1 - s)) ^ 2)
	}
	
	# some sanity bounding
	errors <- pmin(errors, sum(y ^ 2))
	
	errors <- errors / n; # normalize loo-cv error
	
	# determine relevant dimension
	if(dim_rest < 1)
	{
		rd <- which.min(errors[1:(n*dim_rest)])
	}
	else
	{
		rd <- which.min(errors[1:n])
	}
	
	if(est_noise || alldim)
	{
		# if noise should be estimated, estimation of the y's is needed
		est_y <- TRUE
	}
	
	ret <- list(rd = rd, err = errors[1:(n-1)], kpc = z, eigvec = U, eigval = E$values, tcm = FALSE)
	
	if(est_y)
	{
		# find out whether this is a regression or classification problem
		regression <- isregression(y, regression)
		if(!regression)
		{
			# classification problem
			Yh <- sign(Yh)
			Yh[Yh == 0] <- 1
		}
		# additionally return estimated ys
		yh <- Yh[, rd, drop = FALSE]
		ret <- c(ret, list(yh = yh))
		# return also all estimated ys for other dimensions?
		if(alldim)
		{
			ret <- c(ret, list(Yh = Yh))
		}
		# estimate noise?
		if(est_noise)
		{
			# estimate noise
			ret <- c(ret, list(noise = estnoise(y, yh, regression, nmse)))
		}
	}
	
	return(ret)
}

