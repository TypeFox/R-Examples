`selectmodel` <-
function(X, y, kernel = rbfkernel, est_y = FALSE, ydist = FALSE, est_noise = FALSE, regression = FALSE, nmse = TRUE, tcm = TRUE, Xname = "X", ...)
{
	# test if X is a matrix
	if(!is.matrix(X))
	{
		print("X must be a matrix containing samples in its columns")
		return()
	}
	# test if y has same dimension as X has rows
	if(length(y) != nrow(X))
	{
		print("Length of label vector y must be equal to the number of columns of X")
		return()
	}
	# test if Xname is a string
	if(!is.character(Xname))
	{
		print("Xname must be a string specifying the name of the argument for the data matrix X
			in kernel function")
		return()
	}
	# test if there are any parameters for model selection
	dots <- list(...)
	if(length(dots) < 1)
	{
		print("There must be at least one parameter to vary to do model selection")
		return()
	}
	# test if all arguments in ... have a name
	if(length(dots) != length(names(dots)))
	{
		print("All arguments in ... must have a name, namely the same names 
			as in the kernel function")
		return()
	}
	# test if all parameters in ... are lists (of any type)
	for(i in 1:length(dots))
	{
		if(!is.vector(dots[[i]]))
		{
			print("All arguments in ... must be lists of parameter 
				values to use for model selection")
			return()
		}
	}
	
	# denoised ys have to be calculated to estimate noise level
	if(est_noise)
	{
		est_y <- TRUE
	}

	# determine all parameter combinations
	comb <- expand.grid(dots)
	n <- ncol(comb)
	m <- nrow(comb)
	
	# loglik/cv-error matrix
	errs <- matrix(0, m, length(y) - 1)
	
	# vector with relevant dimensions
	rds <- rep(0, m)
	
	if(ydist)
	{
		est_y <- TRUE
		Yd <- matrix(0, m, length(y))
	}
	
	# save rde results for model with smallest error
	besterr <- Inf
	best <- NULL
	rd <- NULL
	err <- NULL
	kpc <- NULL
	eigvec <- NULL
	eigval <- NULL
	noise <- NULL
	yh <- NULL
	
	# walk over all parameter combinations
	for(i in 1:m)
	{
		# create list with current parameters
		# which will be passed to the kernel function
		p <- list()
		for(j in 1:n)
		{
			p[[names(comb)[j]]] <- comb[[i, j]]
		}
		
		# calculate kernel matrix
		lXname <- list()
		lXname[[Xname]] <- X
		K <- do.call(kernel, c(lXname, p))
		
		# relevant dimension estimation
		r <- rde(K, y, tcm = tcm, est_y = est_y, alldim = ydist, est_noise = est_noise, regression = regression, nmse = nmse)
		
		# save relevant dimension
		rds[i] <- r$rd
		
		# save negativ log-likelihood/cross-validation error
		# for current parameter combination
		errs[i, ] <- r$err
		
		if(ydist)
		{
			Yd[i, ] <- sqrt(apply((r$Yh - matrix(y, length(y), length(y))) ^ 2, 2, sum))
		}
		
		if(errs[i, ][r$rd] < besterr)
		{
			besterr <- errs[i, ][r$rd]
			best <- p
			rd <- r$rd
			err <- errs[i, ]
			kpc <- r$kpc
			eigvec <- r$eigvec
			eigval <- r$eigval
			noise <- r$noise
			yh <- r$yh
		}
	}
	
	ret <- list(rd = rd, rds = rds, err = err, best = best, errs = errs, kpc = kpc, eigvec = eigvec, eigval = eigval, params = dots, tcm = tcm, kernel = kernel, Xname = Xname, X = X, regression = isregression(y, regression))
	if(est_y)
	{
		ret <- c(ret, list(yh = yh))
	}
	if(est_noise)
	{
		ret <- c(ret, list(noise = noise))
	}
	if(ydist)
	{
		ret <- c(ret, list(Yd = Yd))
	}
	return(ret)
}

