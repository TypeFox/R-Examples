checkSymmetricPositiveDefinite <- function(x, name="sigma") {
  if (!isSymmetric(x, tol = sqrt(.Machine$double.eps))) {
		stop(sprintf("%s must be a symmetric matrix", name))
	}
	
	if (NROW(x) != NCOL(x)) {
		stop(sprintf("%s must be a square matrix", name))
	}
	
	if (any(diag(x) <= 0)) {
	  stop(sprintf("%s all diagonal elements must be positive", name))
  }
	
	if (det(x) <= 0) {
		stop(sprintf("%s must be positive definite", name))
	}
}

# Uses partly checks as in mvtnorm:::checkmvArgs!
checkTmvArgs <- function(mean, sigma, lower, upper)
{
	if (is.null(lower) || any(is.na(lower))) 
		stop(sQuote("lower"), " not specified or contains NA")
	if (is.null(upper) || any(is.na(upper))) 
		stop(sQuote("upper"), " not specified or contains NA")
	if (!is.numeric(mean) || !is.vector(mean)) 
		stop(sQuote("mean"), " is not a numeric vector")
	if (is.null(sigma) || any(is.na(sigma))) 
		stop(sQuote("sigma"), " not specified or contains NA")
	
	if (!is.matrix(sigma)) {
		sigma <- as.matrix(sigma)
	}
	
	if (NCOL(lower) != NCOL(upper)) {
		stop("lower and upper have non-conforming size")
	}
	
	checkSymmetricPositiveDefinite(sigma)
	
	if (length(mean) != NROW(sigma)) {
		stop("mean and sigma have non-conforming size")
	}
	
	if (length(lower) != length(mean) || length(upper) != length(mean)) {
		stop("mean, lower and upper must have the same length")
	}
	
	if (any(lower>=upper)) {
		stop("lower must be smaller than or equal to upper (lower<=upper)")
	}
	
	# checked arguments
	cargs <- list(mean=mean, sigma=sigma, lower=lower, upper=upper)
	return(cargs)
}