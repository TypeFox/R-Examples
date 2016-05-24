#' Simple wrapper for Reduce, proceeding from the right.
#'
#' Wrapper for \code{\link[base]{Reduce}} with \code{right} set to \code{TRUE}.
#'
#' @param ... Arguments to be passed to \code{\link[base]{Reduce}}.
#' @family function operators
#' @seealso \code{\link[base]{Reduce}} for code and documentation.
#' @export
Reduce_Right <- function(...) {
  Reduce(..., right = TRUE)
}

