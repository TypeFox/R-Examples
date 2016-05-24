`rbfkernel` <-
function(X, sigma = 1, Y = NULL)
{
	# test if X is a matrix
	if(!is.matrix(X))
	{
		print("X must be a matrix containing samples in its rows")
		return()
	}
	# test if sigma is a number and > 0
	if(length(sigma) != 1 || sigma <= 0)
	{
		print("sigma must be a number > 0 specifying the rbf-kernel width")
		return()
	}
	if(!is.null(Y))
	{
		# test if Y is a matrix
		if(!is.matrix(Y))
		{
			print("Y must be a matrix containing samples in its rows or NULL if it should not be used")
			return()
		}
		# test if vectors in X and Y have same dimension
		if(ncol(X) != ncol(Y))
		{
			print("The samples in the rows of X and Y must be of same dimension")
			return()
		}
	}
	
	n <- nrow(X) # number of samples in X
	
	if(is.null(Y))
	{
		# calculate distance matrix
		XtX <- tcrossprod(X)
		XX <- matrix(1, n) %*% diag(XtX)
		D <- XX - 2*XtX + t(XX) # distance matrix
	}
	else
	{
		m <- nrow(Y) # number of samples in Y
		# calculate distance matrix (between vectors of X and Y)
		XX <- matrix(apply(X ^ 2, 1, sum), n, m)
		YY <- matrix(apply(Y ^ 2, 1, sum), n, m, byrow = TRUE)
		XY <- tcrossprod(X, Y)
		D <- XX - 2*XY + YY
	}
	
	# calculate rbf-kernel matrix
	K <- exp(-D/(2*sigma))
	
	return(K)
}

