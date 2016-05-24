#' Removes module from cache, applying same logic as \code{\link{find.first.R}} to find and remove it
#'
#' @param module name of a module, same as the one used in the \code{\link{lrequire}} method, that will be removed
#'               from the cache, such that the next time the \code{module} is \code{\link{lrequire}}'d, it will be
#'               read and executed.
#'
#' @param character.only a logical value, defaulted to FALSE, that permits an unquoted name to be \code{lrequire}-d.
#'                       Set this to TRUE when passing a variable to \code{lrequire}, requiring a quoted string.
#'
#' @return boolean value yielding success of removal from the cache
#' @export
#'
#' @examples
#' remove.from.module.cache(variables)
remove.from.module.cache <- function(module, character.only = FALSE) {
  # allow module to be specified without quotes
  if (character.only || length(substitute(module)) > 1) { # an expression was passed, just evaluate it
    module <- as.character(module)
  } else {
    module <- as.character(substitute(module))
  }

  filename <- find.first.R(module, character.only = TRUE) # at this point, `module' is converted to a character
  if (!is.na(filename)) {
    filename.mtime <- paste0(filename, '.mtime')
    remove(list=c(filename, filename.mtime), envir=module.cache)
  }

  return (!is.na(filename))
}
