###########################
## Compute discrepancy
###########################

dissvar <- function(diss, weights=NULL, squared=FALSE) {
	if (squared) {
		diss <- diss^2
	}
	if (is.null(weights)) {
		if (inherits(diss, "dist")) {
			return(sum(diss)/(attr(diss, "Size")^2))
		} else if (is.matrix(diss)) {
			return(sum(diss)/(2*(nrow(diss)^2)))
		} else {
			stop("diss argument should be a dist object or a dissimilarity matrix")
		}
	} 
	else {
		isdist <- inherits(diss, "dist")
		if (isdist) {
			n <- attr(diss, "Size")
		} else if (is.matrix(diss)) {
			n <- nrow(diss)
		} else {
			stop("diss argument should be a dist object or a dissimilarity matrix")
		}
		if(is.null(weights)) {
			weights <- rep(1, n)
		}
		dvar <- .Call(TMR_tmrWeightedInertiaDist, diss, as.integer(n), as.integer(isdist), as.integer(1:n), as.double(weights), as.integer(TRUE))
		return(dvar)
	}
}