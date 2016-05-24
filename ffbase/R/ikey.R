#' Creates a unique integer key for unique combinations of rows of an \code{ffdf}
#'
#' Creates a unique integer key for unique combinations of rows of an \code{ffdf}. In database terms this would correspond to a primary or foreign key.\cr
#' Orders the \code{ffdf} decreasingly alongside the columns with NA's as last in the order and creates the integer key. \cr
#'
#' @export
#' @example ../examples/key.R
#' @param x an \code{ffdf}
#' @param ... other parameters passed on to chunk
#' @return An integer ff vector of the same length as the number of rows in \code{x} with unique values for each unique row
ikey <- function(x, ...){
	o <- ffdforder(x, na.last = TRUE, decreasing = FALSE)
	xo <- x[o, , drop=TRUE]
	idx <- duplicated(xo, ...)
	idx <- ffwhich(idx, idx == FALSE, ...)
	res <- ff(0, vmode="integer", length = nrow(x))
	res[idx] <- ff(1, vmode="integer", length = length(idx))
	res <- cumsum(res)
	res[fforder(o)]
}


