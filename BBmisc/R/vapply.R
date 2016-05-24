#' Apply a function with a predefined return value
#'
#' @description
#' These are just wrappers around \code{\link[base]{vapply}} with
#' argument \code{FUN.VALUE} set.
#' The function is expected to return a single \code{logical}, \code{integer},
#' \code{numeric} or \code{character} value, depending on the second letter
#' of the function name.
#'
#' @param x [\code{vector} or \code{list}]\cr
#'   Object to apply function on.
#' @param fun [\code{function}]\cr
#'   Function to apply on each element of \code{x}.
#' @param ... [\code{ANY}]\cr
#'   Additional arguments for \code{fun}.
#' @param use.names [\code{logical(1)}]\cr
#'   Should result be named?
#'   Default is \code{TRUE}.
#' @export
vlapply = function(x, fun, ..., use.names = TRUE) {
  vapply(X = x, FUN = fun, ..., FUN.VALUE = NA, USE.NAMES = use.names)
}

#' @rdname vlapply
#' @export
viapply = function(x, fun, ..., use.names = TRUE) {
  vapply(X = x, FUN = fun, ..., FUN.VALUE = NA_integer_, USE.NAMES = use.names)
}

#' @rdname vlapply
#' @export
vnapply = function(x, fun, ..., use.names = TRUE) {
  vapply(X = x, FUN = fun, ..., FUN.VALUE = NA_real_, USE.NAMES = use.names)
}

#' @rdname vlapply
#' @export
vcapply = function(x, fun, ..., use.names = TRUE) {
  vapply(X = x, FUN = fun, ..., FUN.VALUE = NA_character_, USE.NAMES = use.names)
}
