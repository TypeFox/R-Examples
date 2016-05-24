#' @method [ lfactor
#' @export
"[.lfactor" <- function(x,i, ..., drop=FALSE) {
  xo <- x
  class(xo) <- "factor"
  if(!missing(i)) {
    res <- xo[i,...,drop=FALSE]
  } else {
    res <- xo[...,drop=FALSE]
  }
  attr(res,"llevels") <- llevels(x)
  class(res) <- c("lfactor", "factor")
  if(drop) {
    res <- droplevels(res)
  }
  res
}
