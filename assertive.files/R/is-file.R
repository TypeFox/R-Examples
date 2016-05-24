#' Is the path a directory?
#' Checks to see if the input path is a directory.
#' 
#' @param x File paths.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @param .xname Not intended to be used directly.
#' @return \code{is_dir} returns \code{TRUE} if and only if the input 
#' path is a directory that exists, as determined by \code{file.info}.
#' @examples
#' assert_all_are_dirs(R.home())
#' @importFrom assertive.base is_true
#' @export
is_dir <- function(x, .xname = get_name_in_parent(x))
{  
  x <- coerce_to(x, "character", .xname)
  call_and_name(
    function(x) 
    {
      ok <- file.info(x)$isdir
      causes <- ifelse(
        is.na(ok),
        "nonexistent",
        ifelse(ok, "", "file")
      )
      ok <- is_true(ok) 
      set_cause(ok, causes)
    }, 
    x
  )
}

#' Does the file exist?
#'
#' Checks to see if the input files exist.
#'
#' @param x Input to check.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @param .xname Not intended to be used directly.
#' @return \code{is_existing_file} wraps \code{file.exists}, showing
#' the names of the inputs in the answer.   \code{assert_*_are_existing_files} 
#' return nothing but throws an error if \code{is_existing_file} returns
#' \code{FALSE}.
#' @note Trailing slashes are removed from paths to avoid a lot of false 
#' negatives by the underlying function \code{file.exists}.
#' @seealso \code{\link[base]{file.exists}}.
#' @examples
#' assert_all_are_existing_files(dir())
#' # These examples should fail.
#' assertive.base::dont_stop(
#'   assert_all_are_existing_files("not an existing file (probably)")
#' )
#' @export
is_existing_file <- function(x, .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "character", .xname)
  # file.exists returns FALSE under Windows when there is a trailing slash
  x <- sub("[\\/]+$", "", x)
  call_and_name(
    function(x)
    {
      ok <- file.exists(x)
      set_cause(ok, ifelse(ok, "", "nonexistent"))
    }, 
    x
  )
}

#' Is the file accessible?
#'
#' Checks to see if the input files can be executed/read/written to.
#'
#' @param x Input to check.
#' @param warn_about_windows Logical.  If \code{TRUE}, then calling the function
#' under Windows will throw a warning about the problems with 
#' \code{\link[base]{file.access}}.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @param .xname Not intended to be used directly.
#' @return \code{is_executable_file} wraps \code{file.access}, showing
#' the names of the inputs in the answer.   \code{assert_is_executable_file} 
#' returns nothing but throws an error if \code{is_executable_file} returns
#' \code{FALSE}.
#' @seealso \code{\link[base]{file.access}}.
#' @examples
#' files <- dir()
#' is_readable_file(files)
#' is_writable_file(files, warn_about_windows = FALSE)
#' is_executable_file(files, warn_about_windows = FALSE)
#' @export
is_executable_file <- function(x, warn_about_windows = TRUE, .xname = get_name_in_parent(x))
{
  if(warn_about_windows)
  {
    warn_about_file.access_under_windows()
  }
  x <- coerce_to(x, "character", .xname)
  call_and_name(
    function(x)
    {
      ok <- file.access(x, mode = 1) == 0L
      set_cause(
        ok, 
        ifelse(file.exists(x), "unexecutable", "nonexistent")
      )
    }, 
    x
  )
}

#' @rdname is_executable_file
#' @export
is_ex_file <- function(x)
{
  .Defunct("is_executable_file")
  is_executable_file(x)
}

#' Is the directory a known R library?
#' 
#' Checks to see if the input directories are known R libraries.
#' 
#' @param x Directory paths
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @param .xname Not intended to be used directly.
#' @note Input paths are converted to character, and then normalized using
#' \code{normalizePaths}.
#' @return \code{is_library} returns \code{TRUE} if and only if the input
#' paths are known R package libraries.  That is, they must be paths
#' returned by \code{.libPaths}.
#' @examples
#' is_library(c(.libPaths(), R.home()))
#' @export
is_library <- function(x, .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "character", .xname)
  call_and_name(
    function(x) 
    {
      std_x <- normalizePath(path.expand(x), winslash = "/", mustWork = FALSE)
      set_cause(std_x %in% .libPaths(), "not a lib")
    }, 
    x
  )
}

#' @rdname is_executable_file
#' @export
is_readable_file <- function(x, warn_about_windows = TRUE, .xname = get_name_in_parent(x))
{
  if(warn_about_windows)
  {
    warn_about_file.access_under_windows()
  }
  x <- coerce_to(x, "character", get_name_in_parent(x))  
  call_and_name(
    function(x)
    {
      ok <- file.access(x, mode = 4) == 0L
      set_cause(
        ok, 
        ifelse(file.exists(x), "unreadable", "nonexistent")
      )
    }, 
    x
  )
}

#' @rdname is_executable_file
#' @export
is_writable_file <- function(x, warn_about_windows = TRUE, .xname = get_name_in_parent(x))
{
  if(warn_about_windows)
  {
    warn_about_file.access_under_windows()
  }
  x <- coerce_to(x, "character", .xname)
  call_and_name(
    function(x)
    {
      ok <- file.access(x, mode = 2) == 0L
      set_cause(
        ok, 
        ifelse(file.exists(x, "unwritable", "nonexistent"))
      )
    }, 
    x
  )
}

#' Warn about file.access under Windows
#' 
#' If the OS is Windows, throw a warning about using 
#' \code{\link[base]{file.access}}.
#' @return Nothing. Invoked for the side effect of throwing a warning under 
#' Windows.
#' @seealso \code{\link[base]{file.access}}
#' @examples
#' \dontrun{
#' assertive.base::dont_stop(warn_about_file.access_under_windows())
#' }
#' @noRd
warn_about_file.access_under_windows <- function()
{ 
  if(.Platform$OS.type == "windows") # as per is_windows, without the dependency
  {
    warning(
      "This function depends on file.access, which can give unexpected results under Windows."
    )
  }
}
