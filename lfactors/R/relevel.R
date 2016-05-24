#' @method relevel lfactor
#' @export
relevel.lfactor <- function(x, ref, ...) {
	xf <- x
	class(xf) <- "factor"
	res <- relevel(xf, ref, ...)
	ll <- llevels(x)
	ll <- c(ll[levels(x) == ref], ll[levels(x) != ref])
  attr(res,"llevels") <- ll
  class(res) <- c("lfactor", "factor")
  res
}