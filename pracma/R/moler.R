##
##  m o l e r . R
##


moler <- function(n) {
	if (length(n) != 1 || n != round(n))
		stop("Argument 'n' must be an integer.")
	if (n <= 0) return(c())

	A <- matrix(0, nrow = n, ncol = n)
	for (i in 1:n) {
		A[i, 1:i] <- (1:i) - 2
	}
	A <- A + t(A)
	diag(A) <- 1:n
	A
}
