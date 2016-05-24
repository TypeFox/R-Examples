#' Split a path into components
#'
#' The first normalized path is split on forward and backward slashes and its components returned as
#' character vector. The drive or network home are extracted separately on windows systems and
#' empty on all other systems.
#'
#' @param path [\code{character(1)}]\cr
#'  Path to split as string
#' @return \code{named list}: List with components \dQuote{drive} (\code{character(1)}
#'  and \dQuote{path} (\code{character(n)}.
#' @export
splitPath = function(path) {
  assertString(path)
  path = normalizePath(path, mustWork = FALSE)
  if (isWindows()) {
    pattern = "^([[:alpha:]]:)|(\\\\[[:alnum:]]+)"
    m = regexpr(pattern, path)
    if (length(m) == 1L && m == -1L)
      stop("Error extracting the drive letter")
    drive = regmatches(path, m)
    regmatches(path, m) = ""
  } else {
    drive = character(0L)
  }
  list(drive = drive, path = Filter(nzchar, strsplit(path, "[/\\]+")[[1L]]))
}
