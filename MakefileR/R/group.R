#' Creates a group of Makefile items
#'
#' Helps separating similar rules.
#'
#' Use the
#' \code{\link[base]{c}} function or the \code{\link[base]{+}} operator
#' to append groups to other groups and Makefiles (thus creating nested groups).
#'
#' @param ... \code{[MakefileR]}\cr Items created by \code{\link{make_rule}} or other \code{make_}
#'   functions
#' @param .dots \code{[list]}\cr Further rules in addition to \code{...}
#' @param sep \code{[character(1)]}\cr Separator between group items,
#'   \code{NULL} (the default) means no separator.
#' @return An object of class \code{MakefileR_group}
#' @seealso \code{\link{c.MakefileR_group}}
#' @family items
#'
#' @examples
#' makefile(make_rule("all", c("first_target", "second_target")))
#'
#' @references \url{https://www.gnu.org/software/make/manual/}
#'
#' @export
make_group <- function(..., .dots = NULL, sep = NULL) {
  rules <- c(list(...), .dots)
  if (!all(vapply(rules, inherits, logical(1), "MakefileR"))) {
    stop("All members of the group must inherit from class MakefileR")
  }
  structure(list(rules = rules, sep = sep),
            class = c("MakefileR_group", "MakefileR"))
}

#' @export
format.MakefileR_group <- function(x, ...) {
  lapply(x$rules, format) %>%
    Reduce(
      f = function(y, z) c(y, if (length(y) == 0L) NULL else x$sep, z),
      init = character())
}

#' Concatenation of rules
#'
#' Rules can be appended to groups and Makefiles using the
#' \code{\link[base]{c}} function or the \code{\link[base]{+}} operator.
#'
#' @param ...,x,y \code{[MakefileR]}\cr Rules, the first
#'   (\code{x} or the first element of \code{...})
#'   must be of class \code{MakefileR_group}
#'   (created by \code{\link{make_group}} or \code{\link{makefile}})
#' @param recursive \code{[any]}\cr Unused
#'
#' @rdname Concatenation
#' @export
#' @examples
#' c(make_group(sep = ""),
#'   make_group(make_comment("Dummy targets"),
#'              make_rule(".FORCE"), make_rule(".SILENT")),
#'   make_group(make_comment("Definitions"),
#'              make_def("A", "a")))
#'
#' makefile() + (make_group() + make_comment("Definitions") + make_def("A", "a"))
c.MakefileR_group <- function(..., recursive = FALSE) {
  rules = list(...)
  first_rule <- rules[[1L]]
  other_rules <- rules[-1L]
  structure(
    make_group(.dots = c(first_rule$rules, other_rules), sep = first_rule$sep),
    class = class(first_rule)
  )
}

#' @export
#' @rdname Concatenation
`+.MakefileR_group` <- function(x, y) c(x, y)
