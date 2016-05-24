
# compute A^n where A is a square matrix, allowing non-integer and negative powers
# For integer values, see the technique in...
#http://en.wikipedia.org/wiki/Exponentiation_by_squaring

mpower <-
function(A,n){
	is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
	if(!ncol(A)==nrow(A)) stop("Input must be a square matrix")
	p <- nrow(A)
	if (n==0) return (diag(p))
	if (n==1) return (A)

	if (n < 0 ) {
		A <- solve(A)
		n <- abs(n)
	}
	if (is.wholenumber(n)) {
	  result <- diag(p)
	  while (n > 0) {
	    if (n %% 2 != 0) {
	      result <- result %*% A
	      n <- n - 1
	    }
	    A <- A %*% A
	    n <- n / 2
	  }
	}
	else {
		if( any( abs(A-t(A)) > 100*.Machine$double.eps  ) ) 
      stop("Input must be a symmetric matrix for non-integer powers")
		result <- with(eigen(A), vectors %*% (values^n * t(vectors)))
	}
	dimnames(result) <- dimnames(A)
	return(result)
}

"%^%" <- function(A,n) mpower(A,n)
