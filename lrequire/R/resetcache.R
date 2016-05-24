#' Resets the module cache, ensuring files are loaded on next require
#'
#' Note the chace contains other hidden variables, kept in the cache (or environment). These are not removed,
#' and removing them will conflict with the ability of \code{\link{lrequire}} to perform properly.
#'
#' @return Nothing is returned
#' @export
#'
#' @examples
#' reset.module.cache()
reset.module.cache <- function() {
  remove(list = ls(module.cache), envir = module.cache)
}
