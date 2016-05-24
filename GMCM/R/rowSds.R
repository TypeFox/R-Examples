#' @rdname colSds
#' @return \code{rowSds} returns a numeric vector of length \code{n}.
#' @examples
#' y <- matrix(rnorm(50), 10, 5)
#' GMCM:::rowSds(y)
#' @keywords internal
rowSds <- function(x) {
  ans <- rowSdsArma(x)
  dim(ans) <- NULL
  names(ans) <- rownames(x)
  return(ans)
}

