matrix.bigq <- function(data=NA,nrow=1, ncol=1, byrow=FALSE,
                        dimnames =NULL, den = NA, ...)
{
    if(!is.null(dimnames))
        warning("'dimnames' are not implemented for this class")
    .Call(as_matrixq,
          data, as.integer(nrow), as.integer(ncol),
	  as.integer(byrow), den)
}


as.matrix.bigq <- function(x, ...) {
    if(is.matrix(x) || is.data.frame(x) || length(dim(x)) == 2L) {
	d <- dim(x)
	n <- d[1L]
	p <- d[2L]
    } else {
        n <- length(x)
        p <- 1L
    }
    matrix.bigq(x, nrow=n, ncol=p)
}

as.vector.bigq <- function(x, mode="any") {
    if(mode == "list") {
	## "easy"; TODO: gmpToListQ()  withOUT matrix
	.Call(gmpMatToListQ, matrix.bigq(x, nrow=length(x), ncol=1L), 1L)
    } else if (mode == "any") {
	attr(x,"nrow") <- NULL
	x
    } else as.vector(as.double(x), mode=mode)
}


t.bigq <- function(x) .Call(bigq_transposeR, x)

cbind.bigq <- function(..., recursive = FALSE)
    .Call(bigrational_cbind, list(...))

rbind.bigq <- function(..., recursive = FALSE)
    .Call(bigrational_rbind, list(...))

apply.bigq <- function(X, MARGIN, FUN, ...)
{
  ## change matrix to a list:
  X <- .Call(gmpMatToListQ, X, as.integer(MARGIN))
  ## then use std lapply() and convert back to vector:
  .Call(bigrational_c, lapply(X, FUN, ...))
}


