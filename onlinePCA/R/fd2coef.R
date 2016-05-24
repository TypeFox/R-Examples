fd2coef <- function (x, basis, byrow = TRUE) 
{
 	if (!is.matrix(x))
 		return(as.vector(basis$S %*% x))
    if (byrow) 
        return(tcrossprod(x, basis$S))
    return(basis$S %*% x)
}
