.rowsum <- function(x, group, reorder=TRUE, na.rm = FALSE, ...) {

  x <- makeNd (x, 2)
  old <- attributes (x)

  x <- base:::rowsum (x, group = group, reorder = reorder, na.rm = na.rm)

  old$old [[1]]$dim [1] <- nrow (x)
  old$old [[1]]$dimnames [1] <- list (rownames (x))

  mostattributes (x) <- old

  x <- restoredim (x)

  x
}

.test (.rowsum) <- function (){
  rowsum (v, seq_along (v))
}

##' rowsum for arrays
##'
##' This function extends the base function \code{\link[base]{rowsum}}.
##' @param x array to be \code{rowsum}med
##' @param group grouping variable (integer or factor) indicating groups of samples. \code{}
##' @param reorder should the groups be ordered? see \code{\link[base]{rowsum}}
##' @param na.rm shoud \code{NA}s be removed?
##' @param ... ignored
##' @return like \code{\link[base]{rowsum}}, but further dimensions of the array are preserved.
##' @author Claudia Beleites
##' @seealso \code{\link[base]{rowsum}}
##' @keywords array algebra arith
##' @docType methods
##' @rdname rowsum
##' @export
setMethod ("rowsum", signature = signature (x = "array"), .rowsum)

#TODO: test, examples
