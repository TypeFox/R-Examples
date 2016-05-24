##
##  k r o n . R  Kronecker product
##


kron <- function(a, b) {	
	if (length(a) == 0 || length(b) == 0) return(c())
	if (!(is.numeric(a) || is.complex(a)) ||
		!(is.numeric(b) || is.complex(b)))
		stop("Arguments 'a' and 'b' must be real/complex vectors/matrices.")
	
	kronecker(a, b)
}
