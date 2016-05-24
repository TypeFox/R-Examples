`kpr` <-
function(model, X = NULL, Xname = "X", Yname = "Y", kernel = NULL, regression = TRUE, ...)
{
	# test if X is a matrix
	if(!is.null(X))
	{
		if(!is.matrix(X))
		{
			print("X must be a matrix containing samples in its rows")
			return()
		}
	}
	else
	{
		if(is.null(model$X))
		{
			print("You must specify a sample matrix X")
			return()
		}
		else
		{
			X <- model$X
		}
	}
	# test if a kernel was specified (by user or in model)
	if(is.null(kernel))
	{
		if(is.null(model$kernel))
		{
			print("You must specify a kernel function")
			return()
		}
		else
		{
			kernel <- model$kernel
		}
	}
	# test if Xname is specified in model
	if(!is.null(model$Xname))
	{
		Xname <- model$Xname
	}
	# test if regression is specified in model
	if(!is.null(model$regression))
	{
		regression <- model$regression
	}
	
	n <- length(model$eigval)
	# divide eigenvectors through eigenvalues of kernelmatrix
	alpha <- (model$eigvec[,1:model$rd] / matrix(model$eigval[1:model$rd], n, model$rd, byrow=TRUE)) %*%
			model$kpc[1:model$rd, , drop = FALSE]

	# make copy of X (avoids naming conflict in returned function)
	Xc <- X

	return(
			# prediction function for kernel projection regression
			function(X)
			{
				lXY <- list()
				lXY[[Xname]] <- Xc
				lXY[[Yname]] <- as.matrix(X)
				K <- do.call(kernel, c(lXY, model$best, list(...)))
				yh <- apply(matrix(alpha, n, ncol(K)) * K, 2, sum)
				
				if(!regression)
				{
					# classification problem
					yh <- sign(yh)
					yh[yh == 0] <- 1
				}
				return(yh)
			}
		)
}

