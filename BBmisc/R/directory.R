#' Is one / are several files a directory?
#'
#' If a file does not exist, \code{FALSE} is returned.
#'
#' @param ... [\code{character(1)}]\cr
#'   File names, all strings.
#' @return [\code{logical}].
#' @export
#' @examples
#' print(isDirectory(tempdir()))
#' print(isDirectory(tempfile()))
isDirectory = function(...) {
  paths = c(...)
  if (.Platform$OS.type == "windows" && getRversion() < "3.0.2")
    paths = sub("^([[:alpha:]]:)[/\\]*$", "\\1//", paths)
  x = file.info(paths)$isdir
  !is.na(x) & x
}

#' Is one / are several directories empty?
#'
#' If file does not exist or is not a directory, \code{FALSE} is returned.
#'
#' @param ... [\code{character(1)}]\cr
#'   Directory names, all strings.
#' @return [\code{logical}].
#' @export
#' @examples
#' print(isEmptyDirectory(tempdir()))
#' print(isEmptyDirectory(tempfile()))
isEmptyDirectory = function(...) {
  vapply(list(...), FUN.VALUE = TRUE, FUN = function(x) {
    isDirectory(x) && length(list.files(x, all.files = TRUE, include.dirs = TRUE)) == 2L
  })
}
