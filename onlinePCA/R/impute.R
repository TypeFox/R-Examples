impute <- function(lambda, U, x, center, tol = 1e-7)
{
	if (!missing(center))
		x <- x - center
	na <- which(is.na(x))
    if (length(na) == 0) 
    	return(x)
    if (length(na) == length(x)) 
		stop("x contains only NAs")

 	A <- U %*% diag(sqrt(lambda))
		if (nrow(U) - length(na) >= ncol(U)) {
        	ginvAx.nona <- suppressWarnings(lsfit(A[-na, , drop = FALSE], 
                x[-na], intercept = FALSE)$coefficients)
            x[na] <- A[na, , drop = FALSE] %*% ginvAx.nona
        }
        else {
            svdA <- svd(A)
            pos <- svdA$d > tol
            if (!all(pos)) {
                Ainv <- tcrossprod(svdA$v[, pos, drop = FALSE] %*% 
                  diag(1/svdA$d[pos]), svdA$u[, pos, drop = FALSE])
                x[na] <- A[na, , drop = FALSE] %*% Ainv %*% x[-na]
            }
            else x[na] <- 0
        }
    
	return(x)	
}
