trapzmat <- function(X,Y,delta=1,wt=rep(1,n)) {
#TRAPZMAT integrates the products of two matrices of values
#   using the trapezoidal rule, assuming equal spacing
#  X is the first  matrix of values
#  Y is the second matrix of values
#  DELTA is the spacing between argument values (one by default)
#  WT is a vector of weights (ones by default)
#
#  XtWY is a matrix of integral estimates, number of rows equal to
#  number of col of X, number of cols equal to number of cols of Y

	X <- as.matrix(X)
	Y <- as.matrix(Y)
	
	n <- dim(X)[1]

	if (dim(Y)[1] != n) {
    	stop("X and Y do not have same number of rows.")
	}

	if (length(wt) != n) {
    	stop("X and WT do not have same number of rows.")
	}

	if (delta <= 0) {
    	stop("DELTA is not a positive value.")
	}

	wt[c(1,n)] <- wt[c(1,n)]/2
	wt <- wt*delta

	X <- X*outer(wt,rep(1,dim(X)[2]))
	XtWY <- crossprod(X,Y)
	return(XtWY)
}

