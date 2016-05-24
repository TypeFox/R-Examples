`rde_tcm` <-
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

	# calculate eigen decomposition of K
	E <- eigen(K, symmetric = TRUE)
	U <- E$vectors
	
	y <- matrix(y, n, 1) # y as a column vector
	
	z <- t(U) %*% y # kernel pca coefficients
	z2 <- z ^ 2 # squared kernel pca coefficients
	
	# now calculate negative log-likelihood
	
	# matrices to sum z_i^2 (will be used to calculate the variances)
	sum_matrix_1 <- matrix(0, n, n)
	sum_matrix_2 <- sum_matrix_1
	ltm <- lower.tri(sum_matrix_1, diag = TRUE)
	sum_matrix_1[ltm] <- 1
	sum_matrix_2[!ltm] <- 1
	# drop last line because log-likelihood isn't defined for d = n (d = assumed relevant dimension)
	sum_matrix_1 <- sum_matrix_1[1:(n-1),]
	sum_matrix_2 <- sum_matrix_2[1:(n-1),]
	# vector containing dimension numbers from 1 to (n-1)
	d <- matrix(1:(n-1), n - 1, 1)
	# variances for first part (relevant dimension)
	sigma_1 <- (matrix(1, n - 1, 1) / d) * (sum_matrix_1 %*% z2)
	# variances for second part (noise)
	n_d = n - d
	sigma_2 <- (matrix(1, n - 1, 1) / n_d) * (sum_matrix_2 %*% z2)
	# negative log-likelihood
	loglik = (d/n)*log(sigma_1) + (n_d/n)*log(sigma_2)
	
	# determine relevant dimension
	if(dim_rest < 1)
	{
		rd <- which.min(loglik[1:(n*dim_rest)])
	}
	else
	{
		rd <- which.min(loglik[1:(n-1)])
	}
	
	if(est_noise || alldim)
	{
		# if noise should be estimated, estimation of the y's is needed
		est_y <- TRUE
	}
	
	ret <- list(rd = rd, err = loglik, kpc = z, eigvec = U, eigval = E$values, tcm = FALSE)
	
	if(est_y)
	{
		# find out whether this is a regression or classification problem
		regression <- isregression(y, regression)
		# additionally calculate estimated ys
		if(!alldim)
		{
			yh <- denoiselabels(rd, U, z, regression)
		}
		else
		{
			Yh <- denoiselabels(0, U, z, regression)
			yh <- Yh[, rd, drop = FALSE]
		}
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

