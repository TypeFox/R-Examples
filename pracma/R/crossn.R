##
##  c r o s s n . R  Vector Cross Product
##


crossn <- function(A) {
	if (!is.numeric(A))
		stop("Argument 'A' must be numeric.")

	if (is.vector(A) && length(A) == 2) {
		crossA <- c(A[2], -A[1])
	} else {
		if (is.matrix(A) && nrow(A) >= 2 && ncol(A) == nrow(A) + 1) {
			m <- ncol(A)
			crossA <- numeric(m)
			for (i in 1:m)
				crossA[i] <- (-1)^(i+1) * det(A[, -i])
		} else {
			stop("Matrix 'A' must be of size n x (n+1) with n >= 1.")
		}
	}
	return(crossA)
}
