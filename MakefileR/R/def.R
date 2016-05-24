#' Creates a variable definition in a Makefile
#'
#' A variable definition in a \code{Makefile} consists of a variable name
#' and its defition.  Both are separated by the equality sign \code{=}.
#'
#' No quoting is applied to the definition by this function.
#' Currently, both variable and definition are required to be character values
#' of length one.
#'
#' Use the
#' \code{\link[base]{c}} function or the \code{\link[base]{+}} operator
#' to append definitions to groups and Makefiles.
#'
#' @param variable \code{[character(1)]}\cr Variable name
#' @param definition \code{[character(1)]}\cr Definition for this variable
#' @return An object of class \code{MakefileR_def}
#' @seealso \code{\link{makefile}}, \code{\link{make_group}}
#' @family items
#'
#' @examples
#' make_def("R_USER_LIBRARY", .libPaths()[[1L]])
#' makefile() +
#'   make_def("R_USER_LIBRARY", .libPaths()[[1L]])
#'
#' @references \url{https://www.gnu.org/software/make/manual/}
#'
#' @export
make_def <- function(variable, definition) {
  if (length(variable) != 1) {
    stop("variable must be a character value")
  }
  if (length(definition) != 1) {
    stop("definition must be a character value")
  }
  structure(
    list(
      variable = variable,
      definition = definition
    ),
    class = c("MakefileR_def", "MakefileR"))
}

#' @export
format.MakefileR_def <- function(x, ...) {
  sprintf("%s=%s", x$variable, x$definition)
}
