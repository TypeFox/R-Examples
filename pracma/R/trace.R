##
##  t r a c e . R  Matrix trace
##


Trace <- function(a) {
	if (length(a) <= 1) return(a)
	if ((!is.numeric(a) && !is.complex(a)) || !is.matrix(a))
		stop("Argument 'a' must be a real or complex matrix.")
	if (nrow(a) != ncol(a))
		stop("Matrix 'a' must be square.")

	return(sum(diag(a)))
}
