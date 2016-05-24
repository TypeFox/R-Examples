##
##  i n v . R  matrix inverse
##


inv <- function(a) {
	if (length(a) == 0) return(matrix(0, nrow=0, ncol=0))
	if ((!is.numeric(a) && !is.complex(a)) || !is.matrix(a))
		stop("Argument 'a' must be a numeric or complex matrix.")
	if (nrow(a) != ncol(a))
		stop("Matrix 'a' must be square.")

	e <- try(b <- solve(a), silent=TRUE)
	if (class(e) == "try-error") {
		warning("Matrix appears to be singular.")
		b <- rep(Inf, length(a))
		dim(b) <- dim(a)
	}
	return(b)
}
