#' @title Exports to an environment
#' @description This function is a wrapper around \code{\link{assign}} that
#'   exports the contents of a named list to an environment.  The variable names
#'   in the target environment are constructed from the names of the list items
#'   or taken from a separate argument.
#' @param arg.list list of objects, possibly named.
#' @param arg.names names to use for the items in the target environment. Use
#'   the names of \code{arg.list} by default.
#' @param target.env The target environment.  Use the global environment by
#'   default.
#' @return Invisible \code{NULL}.
#' @seealso \link{export}, \link{assign}
#' @examples
#' export.list(list(newly.created.var=5))
#' newly.created.var
#' rm(newly.created.var)
#' @export
#' @author Roland
#' @references \url{http://stackoverflow.com/a/17484932/946850}
export.list <- function(arg.list, arg.names=names(arg.list),
                        target.env=.GlobalEnv) { # nolint
  stopifnot(length(arg.list) == length(arg.names))
  for (i in seq_along(arg.names)) {
    assign(arg.names[i], arg.list[[i]], target.env)
  }
  invisible(NULL)
}

#' @title Exports to an environment
#' @description This function is a wrapper around \code{\link{export.list}} that
#'   exports variables by their name to another environment.
#' @param ... variables to be exported.
#' @param target.env The target environment.  Use the global environment by
#'   default.
#' @return Invisible \code{NULL}.
#' @seealso \link{export.list}, \link{assign}
#' @examples
#' local({
#'   newly.created.var <- 5
#'   export(newly.created.var)
#' })
#' newly.created.var
#' rm(newly.created.var)
#' @export
#' @author Roland
#' @references \url{http://stackoverflow.com/a/17484932/946850}
export <- function(...,
                   target.env=.GlobalEnv) { # nolint
  arg.list <- list(...)
  arg.names <- sapply(match.call()[-1], deparse)
  export.list(arg.list, arg.names, target.env)
}
