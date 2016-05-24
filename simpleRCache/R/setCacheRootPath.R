#' Set Cache Root Path
#'
#' @param path the cache root directory (DEFAULT: a temporary directory, 
#'   but this can be change to a more permanent location) 
#' 
#' @examples 
#' setCacheRootPath()
#' \dontrun{
#' setCacheRootPath("~/.simpleRCache")
#' }
#' 
#' @export
setCacheRootPath <- function(path=tempdir()) {
  if(!file.exists(path)) {
    dir.create(path)
  }

  options("simpleRCacheRoot" = path)
}
