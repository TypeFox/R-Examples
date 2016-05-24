robreg3S <- function(y, x, dummies=NULL, filter=TRUE, alpha=0.20, K=5, ...){
  if (is.null(dummies)) {
    result <- .robreg3S_noDummies(y, x, filter, alpha, ...)
  } else {
    result <- .robreg3S_withDummies(y, x, dummies, filter, alpha, K, ...)
  }
	class(result) <- "robreg3S"
	return( result )
}



