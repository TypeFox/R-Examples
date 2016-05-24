#' Is a file too big or small?
#' 
#' Checks to see if a file is within a given size range.
#' 
#' @param x Input to check.
#' @param lower Smallest file size allowed, in bytes.
#' @param upper Largest file size allowed, in bytes.
#' @param lower_is_strict If \code{TRUE}, the lower bound is open (strict) 
#' otherwise it is closed.
#' @param upper_is_strict If \code{TRUE}, the upper bound is open (strict)
#' otherwise it is closed.
#' @param na_ignore A logical value.  If \code{FALSE}, \code{NA} values
#' cause an error; otherwise they do not.  Like \code{na.rm} in many
#' stats package functions, except that the position of the failing
#' values does not change.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @param .xname Not intended to be used directly.
#' @return \code{is_empty_file} wraps \code{file.info}, retuning \code{TRUE} 
#' when the input is a file that exists with size zero.  
#'\code{assert_*_are_empty_files} return nothing but throws an error if 
#'\code{is_empty_file} returns \code{FALSE}.
#' @seealso \code{\link[base]{file.info}}.
#' @examples
#' tf <- tempfile()
#' file.create(tf)
#' is_empty_file(tf)
#' cat("some stuff", file = tf)
#' is_non_empty_file(tf)
#' assertive.base::dont_stop(assert_all_file_sizes_are_in_range(tf, lower = 100))
#' unlink(tf)  
#' @importFrom assertive.base is_true
#' @export
is_empty_file <- function(x, .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "character", .xname)
  call_and_name(
    function(x) 
    {
      f_info <- file.info(x)
      ok <- is_true(f_info$size == 0 & !f_info$isdir)
      causes <- ifelse(
        is.na(f_info$size),
        "nonexistent",
        ifelse(
          f_info$size > 0,
          "nonempty", 
          ifelse(f_info$isdir, "dir", "")
        )
      )
      set_cause(ok, causes)
    }, 
    x
  )
}

#' @rdname is_empty_file
#' @export
is_non_empty_file <- function(x, .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "character", .xname)
  call_and_name(
    function(x) 
    {
      f_info <- file.info(x)
      ok <- is_true(f_info$size > 0 & !f_info$isdir)
      causes <- ifelse(
        is.na(f_info$size),
        "nonexistent",
        ifelse(
          f_info$size == 0,
          "empty", 
          ifelse(f_info$isdir, "dir", "")
        )
      )
      set_cause(ok, causes)
    }, 
    x
  )
}

#' @rdname is_empty_file
#' @importFrom assertive.numbers is_in_range
#' @export
is_file_size_in_range <- function(x, lower = 0, upper = Inf, lower_is_strict = FALSE, upper_is_strict = FALSE, .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "character", .xname)
  call_and_name(
    function(x) 
    {
      f_info <- file.info(x)
      ok <- is_in_range(
        f_info$size, 
        lower, 
        upper, 
        lower_is_strict, 
        upper_is_strict
      )
      causes <- cause(ok)
      causes <- gsub("high", "big", causes)
      causes <- gsub("low", "small", causes)
      ok[is.na(f_info$size)] <- FALSE
      causes[is.na(f_info$size)] <- "nonexistent"
      ok[f_info$isdir] <- FALSE
      causes[f_info$isdir] <- "dir"
      set_cause(ok, causes)
    }, 
    x
  )
}


