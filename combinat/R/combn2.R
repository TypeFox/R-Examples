"combn2"<-
function(x, n)
{
#   DATE WRITTEN:  14 April 1994           LAST REVISED:  14 April 1994
#   AUTHOR:  Scott D. Chasalow
#
#   DESCRIPTION:
#         Generate all combinations of the elements of x taken two at a time. 
#         If x is missing,  generate all combinations of 1:n taken two
#         at a time (that is,  the indices of x that would give all 
#         combinations of the elements of x if x with length n had been given).
#         Exactly one of arguments "x" and "n" should be given.
#
	if(!missing(x)) {
		if(!missing(n))
			warning(paste("Only one of arguments x and n allowed;", 
				"argument n was ignored"))
		n <- length(x)
	}
	else if(missing(n))
		stop("Arguments \"x\" and \"n\" both missing")
	if(length(n) > 1) {
		warning(paste("Argument n has", length(n), 
			"elements: only the first used"))
		n <- n[1]
	}
	if(n == 0)
		return(NULL)
	rmat <- array(seq(length = n), c(n, n))	# row(matrix(0,n,n))
	cmat <- t(rmat)	# col(matrix(0,n,n))
	lower.t <- rmat > cmat	# lower.tri(matrix(0,n,n))
	i1 <- cmat[lower.t]
	i2 <- rmat[lower.t]
	if(missing(x))
		cbind(i1, i2)
	else cbind(x[i1], x[i2])
}

