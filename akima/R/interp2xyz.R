##' From an akima() result, produce a 3 column matrix "cbind(x, y, z)"
##'
##' @title From interp() Result, Produce 3-column Matrix
##' @param al a \code{\link{list}} as produced from interp().
##' @param data.frame logical indicating if result should be data.frame or matrix.
##' @return a matrix (or data.frame) with three (3) columns, called "x", "y", "z".
##' @author Martin Maechler, Jan.18, 2013
interp2xyz <- function(al, data.frame = FALSE) {
    stopifnot(is.list(al), identical(names(al), c("x","y","z")))
    xy <- expand.grid(x = al[["x"]], y = al[["y"]], KEEP.OUT.ATTRS=FALSE)
    cbind(if(!data.frame) data.matrix(xy) else xy,
          z = as.vector(al[["z"]]))
}

