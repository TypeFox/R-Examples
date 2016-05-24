#' Partial apply a function, filling in some arguments.
#'
#' Wrapper for \code{\link[pryr]{partial}}.
#'
#' @param ... Arguments to be passed to \code{\link[pryr]{partial}}.
#' @family function operators
#' @seealso \code{\link[pryr]{partial}} for code and documentation.
#' @export
Partial <- function(...) {
  return(pryr::partial(...))
}

