
statList <- function(X, factors, FUN, drop=FALSE, ...) {
	if (!is.numeric(as.matrix(X))) stop("Argument X must be numeric")
	if (!missing(factors)) {
    glist <- split(X, factors, drop=drop)
    result <- lapply(glist, FUN, ...)
  	}
  else result <- list(FUN(X, ...))
  result
}

colMeansList <- function(X, factors, drop=FALSE, ...) {
	statList(X, factors, FUN=colMeans, drop=drop, ...)
}

covList <- function(X, factors, drop=FALSE, ...) {
	statList(X, factors, FUN=cov, drop=drop, ...)
}


