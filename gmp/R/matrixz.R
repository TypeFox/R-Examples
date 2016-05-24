matrix <- function(data=NA, nrow=1, ncol=1, byrow=FALSE, ...)
  UseMethod("matrix")

matrix.default <- function(...) base::matrix(...)
##
## looks "better", but results wrongly, e.g. matrix(1:6, 3) :
## matrix.default <- function(data = NA, nrow = 1, ncol = 1, byrow = FALSE,
##                            dimnames = NULL, ...)
##     ## here, rightly, any "..." will give an error from base::matrix
##     base::matrix(data, nrow=nrow, ncol=ncol, byrow=byrow,
##                  dimnames=dimnames, ...)

matrix.bigz <- function(data=NA, nrow=1, ncol=1, byrow=FALSE, dimnames=NULL,
                        mod=NA, ...)
  {
    if(!is.null(dimnames))
        warning("'dimnames' are not implemented for this class")
    .Call(as_matrixz, data,
          as.integer(nrow), as.integer(ncol), as.integer(byrow),
          mod)
  }


as.matrix.bigz <- function(x, ...)
{
    if(is.matrix(x) || is.data.frame(x) || length(dim(x)) == 2L) {
	d <- dim(x)
	n <- d[1L]
	p <- d[2L]
    } else {
        n <- length(x)
        p <- 1L
    }
    matrix.bigz(x, nrow=n, ncol=p)
}

as.vector.bigz <- function(x, mode="any") {
    if(mode == "list") {
	## "easy"; TODO: gmpToListZ()  withOUT matrix
	.Call(gmpMatToListZ, matrix.bigz(x, nrow=length(x), ncol=1L), 1L)
    } else if (mode == "any") {
	attr(x,"nrow") <- NULL
	x
    } else as.vector(as.double(x), mode=mode)
}

t.bigz <- function(x) .Call(bigint_transposeR, x)


##aperm.bigz <- function(a,perm, resize= TRUE)
##  {
##    dims <- dim(a)
##    if (missing(perm))
##      perm <- c(1,2)
##    if(perm[1] > perm[2])
##      ans = .Call("bigint_transposeR",
##        a,
##       PACKAGE="gmp")
##   else
##     ans = a
##    if(!resize)
##      dim(ans) <- dims
##    ans
##}

## MM: the  attr(*, "nrow")  implementation is hack -- not allowing >= 3d arrays
## --==> advertize  .dimZQ(), .nrowZQ(), etc  instead of the "nrow" attr

##' Is 'x' a "bigz" or "bigq" matrix
##' @return \code{ !is.null(attr(x,"nrow")) } but faster
is.matrixZQ <- function(x) .Call(is_matrix_zq, x)

.dimZQ <- function(x) {# return NULL for non-array {as standard R}
    if(is.null(n <- attr(x,"nrow"))) n else c(n, if(n) length(x)/n else 0L)
}

.dimsetZQ <- function(x,value)
{
  attr(x,"nrow") <-
    if(is.null(value)) { ## dim(m) <- NULL as in standard R:
      NULL
    } else {
      stopifnot(value == (v <- as.integer(value)), v >= 0, prod(v) == length(x))
      v[1L]
    }
  x
}
.nrowZQ <- function(x) attr(x,"nrow")
.ncolZQ <- function(x) if(is.null(n <- attr(x,"nrow"))) n else if(n) length(x)/n else 0L

dim.bigz <- .dimZQ
dim.bigq <- .dimZQ
`dim<-.bigz` <- .dimsetZQ
`dim<-.bigq` <- .dimsetZQ

## not usable as S3, as long as  ncol(), nrow() are not S3 generic
## FIXME? remove these ... ncol() / nrow() go via dim() anyway!
nrow.bigz <- .nrowZQ
nrow.bigq <- .nrowZQ
ncol.bigz <- .ncolZQ
ncol.bigq <- .ncolZQ


cbind.bigz <- function(..., recursive = FALSE)
    .Call(biginteger_cbind, list(...))

rbind.bigz <- function(..., recursive = FALSE)
    .Call(biginteger_rbind, list(...))

apply <- function(X, MARGIN, FUN, ...)
    UseMethod("apply")
apply.default <- function(X, MARGIN, FUN, ...)
    base::apply(X, MARGIN, FUN, ...)

apply.bigz <- function(X, MARGIN, FUN, ...)
{
  ## change matrix to a list:
  X <- .Call(gmpMatToListZ, X, as.integer(MARGIN))
  ## then use std lapply() and convert back to vector:
  .Call(biginteger_c, lapply(X, FUN, ...))
}


