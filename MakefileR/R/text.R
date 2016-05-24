#' Creates a custom Makefile entry
#'
#' For anything else not covered by the other \code{make_} functions, such as
#' the \code{export} statement for exporting Makefile variables.
#'
#' Use the
#' \code{\link[base]{c}} function or the \code{\link[base]{+}} operator
#' to append comments to groups and Makefiles.
#'
#' @param ... \code{[character]}\cr Custom text
#' @return An object of class \code{MakefileR_text}
#' @family items
#'
#' @examples
#' make_text("export SOME_VARIABLE")
#'
#' @references \url{https://www.gnu.org/software/make/manual/}
#'
#' @export
make_text <- function(...) {
  text <- c(...)
  if (length(text) == 0L) {
    stop("At least one element is required for a text")
  }
  structure(list(text = text),
            class = c("MakefileR_text", "MakefileR"))
}

#' @export
format.MakefileR_text <- function(x, ...) {
  x$text
}
