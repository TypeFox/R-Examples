TransMatrix.TPmsm <- function(x) {
	mat <- matrix(0, nrow=3, ncol=3)
	rownames(mat) <- colnames(mat) <- x$state.names
	if ( !is.null(x$n.boot) ) lst <- vector(mode="list", length=3)
	else lst <- vector(mode="list", length=1)
	n <- length(x$time)
	for ( i in 1:length(lst) ) {
		mat[1,1] <- x[[i+1]][n,1]
		mat[1,2] <- x[[i+1]][n,2]
		mat[1,3] <- x[[i+1]][n,3]
		mat[2,2] <- x[[i+1]][n,4]
		mat[2,3] <- x[[i+1]][n,5]
		mat[3,3] <- 1
		lst[[i]] <- mat
	}
	return(lst)
}

TransMatrix.TPCmsm <- function(x) {
	mat <- matrix(0, nrow=3, ncol=3)
	rownames(mat) <- colnames(mat) <- x$state.names
	nx <- length(x$x)
	lst <- vector(mode="list", length=nx)
	if ( !is.null(x$n.boot) ) nm <- 3
	else nm <- 1
	nt <- length(x$time)
	for (i in 1:nx) {
		lst[[i]] <- vector(mode="list", length=nm)
		ix <- which(x$covariate == x$x[i])
		for (j in 1:nm) {
			mat[1,1] <- x[[j+1]][nt, ix, 1]
			mat[1,2] <- x[[j+1]][nt, ix, 2]
			mat[1,3] <- x[[j+1]][nt, ix, 3]
			mat[2,2] <- x[[j+1]][nt, ix, 4]
			mat[2,3] <- x[[j+1]][nt, ix, 5]
			mat[3,3] <- 1
			lst[[i]][[j]] <- mat
		}
	}
	return(lst)
}
