fastmerge <- function(DataMat, w = NULL, tol = 1e-04){
	DataMat <- as.matrix(DataMat)
	p <- ncol(DataMat)
	x <- DataMat[,1]
	DataMat <- DataMat[order(x),]
	n <- length(x)
	if(is.null(w)){ w <- rep(1,n) }
	xx <- round({x - mean(x)}/tol)
	nd <- !duplicated(xx)
	ux <- sort(x[nd])
	uxx <- sort(xx[nd])
	nx <- length(ux)
	if (nx == n) {
		ox <- TRUE
		tmp <- cbind(w, DataMat)
	} else {
		ox <- match(xx, uxx)
		tapply1 <- function(X, INDEX, FUN = NULL, ..., simplify = TRUE){
			sapply(X = unname(split(X, INDEX)), FUN = FUN, ..., simplify = simplify, USE.NAMES = FALSE)
		}
		foo <- function(i, D, q) c(sum(q[i]), colMeans(D[i,,drop = FALSE]))
		tmp <- matrix(
				unlist(
				tapply1(seq_len(n), ox, foo, D = DataMat, q = w),
				use.names = FALSE),
				ncol = p+1, byrow = TRUE
			   )
	}
	tmp <- tmp[order(tmp[,2L]),]
	w <- tmp[, 1L]
	DataMat <- tmp[, -1L]
	return(list(DataMat = DataMat, w = w))
}