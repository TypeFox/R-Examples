#' Returns the current file cache
#'
#' @return environment containing the file cache.
#' @export
#'
#' @examples
#' cache <- get.module.cache()
get.module.cache <- function() {
  return(module.cache)
}
