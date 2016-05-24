#' Trace function
#'
#' Matrix algebra
#' @param m a square matrix
#' tr()

tr <- function(m) {
	if (!is.matrix(m) | (dim(m)[1] != dim(m)[2])) 
        stop("m must be a square matrix")
    return(sum(diag(m)))
}

