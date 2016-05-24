`polykernel` <-
function(X, d, Y = NULL)
{
	# test if X is a matrix
	if(!is.matrix(X))
	{
		print("X must be a matrix containing samples in its rows")
		return()
	}
	# test if d is a number and > 0
	if(length(d) != 1 || d <= 0)
	{
		print("d must be a number > 0 specifying the polynomial kernel degree")
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
		
	n <- ncol(X) # number of samples
	
	if(is.null(Y))
	{
		return((tcrossprod(X) + 1) ^ d) # calculate and return polynomial kernel
	}
	else
	{
		return((tcrossprod(X, Y) + 1) ^ d) # calculate and return polynomial kernel
	}
}

