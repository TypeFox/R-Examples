coef2fd <- function(beta, basis, byrow = TRUE)
{
	if (!is.matrix(beta))
		return(as.vector(basis$B %*% (basis$invsqrtM %*% beta)))
    if (byrow) 
        return(tcrossprod(beta %*% basis$invsqrtM, basis$B))
    return(basis$B %*% (basis$invsqrtM %*% beta))
}
