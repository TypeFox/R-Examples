module.readable <- function(module) {
  # returns if the module is readable and is not a directory
  return (file.access(module, mode=4) == 0 && !file.info(module)$isdir)
}

#' Returns the path of the first found instance of \code{module} in \code{module.path}.
#'
#' A symbol maybe passed instead of a string for readability. If an expression is passed, it must return a string value.
#'
#' @param module a string (or symbol) specifying the \code{module} to search for existance and readability in the current directory,
#'             and if it cannot be found, searches for it in the list of directories specified by \code{module.paths}
#'             and then through the set of paths with the module using a .R extension, if it was not originally
#'             specified.
#'
#' @param character.only a logical value, defaulted to FALSE, that permits an unquoted name to be \code{lrequire}-d.
#'                       Set this to TRUE when passing a variable to \code{lrequire}, requiring a quoted string.
#'
#' @param warn.not.found a logical value, defaulted to TRUE, can be set to not display warning messages when
#'                       module is not found.
#'
#' @return A string consisting of the path the module was first found searching through module.paths.
#' @export
#'
#' @examples
#' hide.not.found.warnings()  # don't warn on files not foudn by find.first.R()
#'
#' # Returns the path to the first found module according to module.paths
#' hello_ex.path <- find.first.R(hello_ex)

find.first.R <- function(module, character.only = FALSE, warn.not.found = TRUE) {
  # Check if module exists in the following directories (in the specified order). First check for module in current
  # directory, then module.R in current directory. If neither exist, check module module in the remaining paths, followed
  # by module.R in those paths.

  # allow module to be specified without quotes
  if (character.only || length(substitute(module)) > 1) { # an expression was passed, just evaluate it
    module <- as.character(module)
  } else {
    module <- as.character(substitute(module))
  }

  # TODO: Need to simplify complex paths using .. and other special approachs

  # Convert any \ characters to /
  module <- gsub('\\\\', '/', module)

  module.paths <- get('.:module.paths', envir = module.cache)
  search.for.r <- (toupper(substr(module, nchar(module)-1, nchar(module))) != '.R')
  files <- c()

  if ((substr(module, 1, 1) == '/') || (substr(module, 2, 2) == ':')) {
    files <- append(files, ifelse(search.for.r, paste(module, 'R', sep='.'), module))
  } else {
    files <- append(files, file.path('.', ifelse(search.for.r, paste(module, 'R', sep='.'), module)))
    files <- append(files, file.path(module.paths, module))

    # Check module.R only if .R was not specified
    if (search.for.r) {
      files <- append(files, file.path(module.paths, paste(module, 'R', sep='.')))
    }
  }

  # Condense multiple references to same directory ('./' --> '')
  files <- gsub('(?<!\\.)\\./', '', files, perl=TRUE)

  for (filepath in files) {
    if ((substr(filepath, 1, 1) == '/') || (substr(filepath, 2, 2) == ':')) {
      # do nothing
    } else {
      filepath <- path.expand(file.path(getwd(), filepath))
      filepath <- gsub('/(?<!\\.)[^/]+/\\.\\.', '', filepath, perl=TRUE)
    }
    if (module.readable(filepath)) {
      return (filepath)
    }
  }

  if (missing(warn.not.found)) { # Only take global value is not specified in parameter call
    warn.not.found <- get('.:warn.not.found', envir = module.cache)
  }
  if (warn.not.found) {
    warning('No module named ', module,
            ifelse(search.for.r, paste0(' or ', module, '.R'), ''),
            ' in: ', paste(c('.', module.paths), collapse='/, '), '/.')
  }
  return(NA)
}

#' Globally hide warnings when modules are not found
#'
#' This is only to be called when handling results manually. This will be overridden by a \code{warn.not.found}
#' parameter explicitly set by either \code{\link{lrequire}}, \code{\link{find.first.R}}.
#'
#' @return nothing is returned
#' @export
#'
#' @examples
#' # Ensure warnings are not displayed when lrequire cannot find the module
#' hide.not.found.warnings()
hide.not.found.warnings <- function() {
  assign('.:warn.not.found', FALSE, envir = module.cache)
}

#' Globally show warnings when modules are not found
#'
#' This is the default behaviour.
#'
#' This is only to be called when handling results manually. This will be overridden by a \code{warn.not.found}
#' parameter explicitly set by either \code{\link{lrequire}} or \code{\link{find.first.R}}.
#'
#' @return nothing is returned
#' @export
#'
#' @examples
#' # Ensure warnings are displayed when lrequire cannot find the module
#' show.not.found.warnings()
show.not.found.warnings <- function() {
  assign('.:warn.not.found', TRUE, envir = module.cache)
}

