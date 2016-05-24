`denoiselabels` <-
function(d, eigvec, kpc, regression = TRUE)
{
	# test if d is a number >= 0
	if(length(d) != 1 || d < 0)
	{
		print("d >= 0 must be the dimension to project the labels to or 0 if the labels should be projected 
			to each dimension and a matrix with all projections should be returned")
		return()
	}
	# test if eigvec is a matrix
	if(!is.matrix(eigvec))
	{
		print("eigvec must be a matrix containing the eigenvectors of the kernel matrix")
		return()
	}
	# test if eigvec is a square matrix
	if(nrow(eigvec) != ncol(eigvec))
	{
		print("eigvec must be a square matrix containing the eigenvectors of the kernel matrix")
		return()
	}
	# test if kpc is a column vector and has correct size
	if(!is.matrix(kpc))
	{
		print("kpc must be a column vector containing the kernel pca coefficients")
		return()
	}
	if(nrow(kpc) != nrow(eigvec) || ncol(kpc) != 1)
	{
		print("kpc must be a column vector of same dimension as eigvec containing the kernel pca coefficients")
		return()
	}

	n <- nrow(eigvec)
	
	if(d != 0)
	{
		# denoising only for dimension d
		yh <- eigvec[,1:d] %*% kpc[1:d, , drop = FALSE]
	}
	else
	{
		# denoising for all dimensions
		yh <- matrix(rep(0, n*n), n)
		yh[, 1] <- eigvec[, 1, drop = FALSE] * kpc[1]
		for(i in 2:n)
		{
			yh[, i] <- yh[, i - 1, drop = FALSE] + eigvec[, i, drop = FALSE] * kpc[i]
		}
	}
	
	if(!regression)
	{
		# classification problem
		yh <- sign(yh)
		yh[yh == 0] <- 1
	}
	return(yh)
}

