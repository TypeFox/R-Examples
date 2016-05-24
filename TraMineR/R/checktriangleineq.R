## Check all possible 3 distance and check for triangle inequality consistency

checktriangleineq <- function(mat, warn=TRUE, indices = FALSE, tol = 1e-7) {
	## Take care to get a matrix
	mat <- dist2matrix(mat)
	ind <- .Call(TMR_checktriangleineq, mat, as.integer(nrow(mat)), as.double(tol))
	if(is.null(ind)){
		return(TRUE)
	}
	if (warn) {
		warning("At least the indices [", ind[1], ", ", ind[2], "] does not respect the triangle inequality when going through ", ind[3])
	}
	if (indices) {
		return(ind)
	}
	return(FALSE)
	
}