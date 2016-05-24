## transform a ridge object to PCA space, given by the SVD of X

pca.ridge <- function(x, ...) {
	if (is.null(x$svd.V)) stop("ridge object must contain svd.V")
	
	# transform a variance covariance matrix S by L
	rot <- function(S,L) L %*% S %*% t(L)
	
	result <- x	
	result$coef <- result$coef %*% result$svd.V
	result$cov <- lapply(result$cov, rot,  t(result$svd.V))
	class(result) <- c("pcaridge", class(x))
	result
}
