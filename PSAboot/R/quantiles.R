#' Return the 25th percentile.
#' 
#' @param x numeric vector.
#' @param na.rm	logical; if true, any NA and NaN's are removed from x before the
#'        quantiles are computed
#' @param ... other parameters passed to \code{\link{quantile}}.
#' @export
#' @return the 25th percetile.
q25 <- function(x, na.rm=FALSE, ...) {
	return(unname(quantile(x, na.rm=na.rm, ...)[2]))
}

#' Returns the 75th percentile.
#' 
#' @param x numeric vector.
#' @param na.rm	logical; if true, any NA and NaN's are removed from x before the
#'        quantiles are computed
#' @param ... other parameters passed to \code{\link{quantile}}.
#' @export
#' @return the 75th percentile.
q75 <- function(x, na.rm=FALSE, ...) {
	return(unname(quantile(x, na.rm=na.rm, ...)[4]))
}
