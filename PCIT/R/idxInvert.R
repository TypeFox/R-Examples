idxInvert <- function(m, idx) {
	if (class(m)=="numeric" | class(m)=="integer") {
		nNodes <- m
	} else {
		nNodes <- try(nrow(m), silent=TRUE)
		if (class(nNodes) == "try-error" | is.null(nNodes)) {
			cat("ERROR: argument 'm' must be a numeric OR an object on which nrow() can be performed.\n\n", geterrmessage())
			return(FALSE)
		}
	}
	
	return(setdiff(1:(nNodes^2), idx))
}
