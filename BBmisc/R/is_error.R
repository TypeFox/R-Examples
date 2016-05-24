#' Is return value of try an exception?
#'
#' Checks if an object is of class \dQuote{try-error} or
#' \dQuote{error}.
#'
#' @param x [any]\cr
#'   Any object, usually the return value of \code{\link[base]{try}},
#' \code{\link[base]{tryCatch}}, or a function which may return a
#' \code{\link[base]{simpleError}}.
#' @return [\code{logical(1)}].
#' @export
#' @examples
#' x = try(stop("foo"))
#' print(is.error(x))
#' x = 1
#' print(is.error(x))
is.error = function(x) {
  inherits(x, c("try-error", "error"))
}
