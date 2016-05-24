#' Is the path on the OS path?
#'
#' Is the specified path on the operating system search path?
#' 
#' @param x An path to check.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @param .xname Not intended to be used directly.
#' @return \code{TRUE} if the specified paths are on the OS search path.
#' @note The OS search path is determined with \code{Sys.getenv("path")}.  For
#' files, the path of the containing folder is checked rather than the path of 
#' the file itself.
#' @examples
#' is_on_os_path(
#'   c(R.home("bin"), R.home("etc"), "a nonexistent path")
#' ) # probably c(TRUE, FALSE, FALSE)
#' @importFrom assertive.base is_false
#' @importFrom assertive.base call_and_name
#' @importFrom assertive.base coerce_to
#' @importFrom assertive.base set_cause
#' @export
is_on_os_path <- function(x, .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "character", .xname)
  sep <- if(is_windows()) ";" else ":"
  # For files, check the containing directory
  # is_dir treats missing directories as FALSE, which we don't want here
  # so just as easy to use underlying file.info
  x <- ifelse(is_false(file.info(x)$isdir), dirname(x), x)
  call_and_name(
    function(x) 
    {
      x <- normalizePath(path.expand(x), mustWork = FALSE)
      paths <- normalizePath(
        strsplit(Sys.getenv("PATH"), sep)[[1]], 
        mustWork = FALSE
      )
      ok <- x %in% paths
      set_cause(ok, ifelse(file.exists(x), "not on path", "nonexistent"))
    }, 
    x
  )  
}
