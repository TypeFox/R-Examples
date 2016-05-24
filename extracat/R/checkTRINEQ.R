## Check all possible 3 distance and check for triangle inequality consistency

trinq <- function(mat, warn=TRUE, indices = FALSE, tol = 1e-7) {
	## Take care to get a matrix
	mat <- as.matrix(mat)
	n <- nrow(mat)
	## If less than 3, no need for a check
	if (n < 3) return(TRUE)
	
	#Checking only upper triangle of the matrix (without diagonal)
	for (i in 1:(n-1)) {
		## upper triangle, starting at i+1
		for (j in (i+1):n) {
			#Distance between i and j
			d <- mat[i, j]
			## Try to find a z point that break triangle inequality
			for (z in 1:n) {
				## Triangle inequality check
				if (d-(mat[i, z] + mat[z, j]) >= tol) {
					## if warn, shout for problem
					if (warn) {
						warning("At least the indices [", i, ", ", j, "] does not respect the triangle inequality when going through ", z)
					}
					if (indices) return(c(i, j, z))
					return(FALSE)
				}
			}
		}
	}
	## No inconstancy found
	return(TRUE)
}