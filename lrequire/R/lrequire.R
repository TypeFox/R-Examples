#' Sources an R module with optional caching for subsequent attempts, exporting specified values
#'
#' \code{lrequire} looks in the current path, and then through a list of predefined paths
#' to search for the given module to source into the current environment, but only making visible
#' specific variables that are "exported" as a list, in a fashion similar
#' to \href{https://nodejs.org/}{node.js}. The caching behaviour can be either suspended or it can
#' re-source files that have changed since the last time the module was cached.
#'
#' @param module a string (or expression) that specifies a module to load, with or without an optional .R extension. If
#'        the module does not exist in the current directory, it searches for the module in directories listed in
#'        \code{module.paths}, first seaching all directories for the named module,
#'        then the module with a .R extension.
#'
#'        \itemize{
#'          \item{./R_modules/}
#'          \item{./lib/}
#'          \item{../R_modules/}
#'          \item{../lib/}
#'          \item{~/.R_modules}
#'        }
#'
#'        All variables exposed in the module will be hidden in the calling environment, except for
#'        what is exposed through module.exports or the exports list variable.
#'
#' @param force.reload a logical value, defaulted to FALSE, that can be set to TRUE to disable caching behavior for
#'                     the module. If the module has already been loaded and cached, setting \code{force.reload} to
#'                     TRUE will re-source the module. Setting it again to FALSE will re-source the module if the
#'                     previous state was TRUE.
#'
#' @param character.only a logical value, defaulted to FALSE, that permits an unquoted name to be \code{lrequire}-d.
#'                       Set this to TRUE when passing a variable to \code{lrequire}, requiring a quoted string.
#'
#' @param warn.not.found a logical value, defaulted to TRUE, can be set to not display warning messages when
#'                       module is not found.
#'
#' @details
#' \code{lrequire} operates in a similar principle to modules in \href{https://nodejs.org/}{node.js} - keeping
#' any variables created in the source module isolated from the calling environment, while exposing a select set
#' of values/parameters. The specific values are exposed by setting a named list element in the \code{exports} variable
#' to the desired value or by assigning \code{module.exports} a value.
#'
#' Note this list exposed in \code{module.exports} should have named items so they can easily be accessed in
#' the calling environment, however that is not necessary if only a single value is being returned.
#'
#' If values are assigned to both \code{module.exports} and \code{exports}, only the values in \code{module.exports}
#' will be exposed to the caller.
#'
#' Caching a long-running operation, such as static data retrieval from a database is a good use of the
#' caching capability of \code{lrequire} during development when the same module is sourced multiple times.
#'
#' During development, files can be reloaded, even if being cached, if they have been modified after the time they
#' were cached. To enable this behaviour, set the variable \code{module.change_code} to 1.
#'
#' To quickly clear lrequire's package environment, unload the package. In RStudio, this can be done by unchecking
#' \code{lrequire} on the Packages tab. You can also execute the following at the R prompt:
#' \code{
#'     detach("package:lrequire", unload=TRUE)
#'     }
#' The next call to \code{library(lrequire)} will ensure it starts off with a clean slate.
#'
#' @return Any values that exist in \code{module.exports} or, if that does not exist, then the
#'         \emph{list} \code{exports}.
#'
#'         If no module is found, \code{NA} is returned.
#'
#' @author Rick Wargo, \email{lrequire@rickwargo.com}
#'
#' @examples
#' hide.not.found.warnings()  # don't warn on files not found by lrequire()
#'
#' # If the module name is in a character vector, use:
#' my.module <- 'myplot'
#' mm <- lrequire(my.module, character.only = TRUE)
#'
#' say.hello.to <- lrequire(hello_ex)
#' # say.hello.to('Rick')  # use the say.hello.to() function that was returned by lrequire()
#'
#' @export
lrequire <- function(module, force.reload = FALSE, character.only = FALSE, warn.not.found = TRUE) {
  # Look for module or module.R in paths and source locally, exposing list associated with module.exports or exports
  # Return NA otherwise
  #
  # TODO (med): Keep an internal list of dependencies and if a module changes that is dynamically reloaded, reload all dependencies
  # TODO (low): Check if module.change_code is set and if so, have a dynamic watcher to reload the module if it changes
  # TODO (low): allow same filename to be cached multiple ways (maybe cache on getpath)

  # allow module to be specified without quotes
  if (character.only || length(substitute(module)) > 1) { # an expression was passed, just evaluate it
    module <- as.character(module)
  } else {
    module <- as.character(substitute(module))
  }

  if (missing(warn.not.found)) { # Only take global value is not specified in parameter call
    warn.not.found <- get('.:warn.not.found', envir = module.cache)
  }

  (function(module, force.reload) {
    filename <- find.first.R(module, character.only = TRUE, warn.not.found = warn.not.found) # at this point, `module' is converted to a character
    if (!is.na(filename)) {
      mtime = file.info(filename)$mtime
      filename.mtime <- paste0(filename, '.mtime')

      if (exists('.:module.change_code', envir=module.cache)) {
        module.change_code <- get('.:module.change_code', envir=module.cache)
      } else {
        module.change_code <- 0
      }

      if (exists(filename.mtime, envir=module.cache) && (module.change_code > 0)) {
        # if change_code is set and the modification time has changed, force the reload
        if (get(filename.mtime, envir=module.cache) != mtime) {
          force.reload <- TRUE
        }
      }

      if (force.reload || !exists(filename, envir=module.cache)) {
        # forcing a reload of the module has never been read
        source(filename, local = TRUE)
        if (exists('module.change_code')) {
          assign('.:module.change_code', module.change_code, envir=module.cache)
        } else {
          module.change_code <- 0
        }

        # Get the value of module.exports or exports and save to later return to the caller
        return.val <- NULL
        if (exists('module.exports') && !is.null(module.exports)) {
          return.val <- module.exports
        } else if (exists('exports') && !is.null(exports)) {
          return.val <- exports
        }

        # Anytime we load a module, we want to save the module and mod time for later checks
        assign(filename, return.val, envir=module.cache)
        assign(filename.mtime, mtime, envir=module.cache)

        return(return.val)
      } else {
        return(get(filename, envir=module.cache))
      }
    }
    return(NA)
  })(module, force.reload)
}
