batchpca <- function(x, q, center, type = c("data","covariance"), 
	byrow = FALSE)
{
	type <- match.arg(type)
	if (type == "data") {
		if (!missing(center))
			x <- x - matrix(center, nrow(x), ncol(x), byrow)
		nu <- ifelse(byrow,0,q)
		nv <- ifelse(byrow,q,0)
		n <- ifelse(byrow,nrow(x),ncol(x))
		svdx <- svds(x/sqrt(n),q,nu,nv)
		return(list(values = svdx$d[1:q]^2, 
				vectors = if (byrow) svdx$v else svdx$u))
	} else if (type == "covariance") {
		eigx <- eigs_sym(x, q)
		return(list(values = eigx$values, vectors = eigx$vectors))
	} else return(NULL)
}
