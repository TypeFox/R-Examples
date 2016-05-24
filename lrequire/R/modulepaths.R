#' Get existing collection of search paths of where to look for modules.
#'
#' @return Vector of paths (strings) that specify folders where to search for module files, in order.
#' @export
#'
#' @examples
#' # Returns a copy of the current module.paths
#' paths <- get.module.paths()
get.module.paths <- function() {
  return (get('.:module.paths', envir = module.cache))
}

#' Append/Insert a path into module.paths, similar to append()
#'
#' Beware, little error checking is done to see if the indexes are valid.
#' Note the item indexes are 1-based.
#'
#' @param value (relative) path to insert into \code{module.paths}
#' @param after location in \code{module.paths} to append, 0 for the beginning, defaults
#'              to -1 which appends to the end
#'
#' @return Nothing is returned.
#' @export
#'
#' @examples
#' # Inserts `../R/lib' as the second path to search for modules
#' append.module.paths('../R/lib', after = 1)
append.module.paths <- function(value, after = -1) {
  module.paths <-get('.:module.paths', envir = module.cache)

  if (after < 0) {  # append to end of path list if no after arg specified
    after = length(module.paths)
  }
  module.paths <- append(module.paths, value, after)
  assign('.:module.paths', module.paths, envir = module.cache)
}

#' Remove one or more paths from \code{module.paths}
#'
#' Beware, little error checking is done to see if the indexes are valid.
#' Note the item indexes are 1-based.
#'
#' @param index index of item to be removed from path
#' @param ...   optional indexes of items to be removed from path
#'
#' @return Nothing is returned.
#' @export
#'
#' @examples
#' # Removes the2nd and 4th items from module.paths
#' remove.module.paths(2, 4)
remove.module.paths <- function(index, ...) {
  indicies <- c(index, c(...))

  module.paths <-get('.:module.paths', envir = module.cache)
  module.paths <- module.paths[-indicies]
  assign('.:module.paths', module.paths, envir = module.cache)
}
