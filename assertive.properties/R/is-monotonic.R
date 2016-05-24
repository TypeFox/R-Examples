#' Is the vector monotonically increasing or decreasing?
#' 
#' Checks to see if the input is monotonically increasing or decreasing.
#' @param x Input to check.
#' @param strictly Logical. If \code{TRUE}, the input is checked for being 
#' strictly monotonic; that is, consecutive values cannot be equal.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @examples 
#' x <- c(1, 2, 2, 1, 3, 2)
#' is_monotonic_increasing(x)
#' is_monotonic_increasing(x, TRUE)
#' is_monotonic_decreasing(x)
#' is_monotonic_decreasing(x, TRUE)
#' 
#' # Also works with, e.g., dates & times
#' is_monotonic_increasing(Sys.time() + x)
#' 
#' # These checks should fail
#' assertive.base::dont_stop({
#'   assert_is_monotonic_increasing(x)
#'   assert_is_monotonic_decreasing(x)
#' })
#' @aliases is_monotonic
#' @export
is_monotonic_increasing <- function(x, strictly = FALSE, .xname = get_name_in_parent(x))
{
  bad_diff_x <- which(
    if(strictly)
    {
      diff(x) <= 0
    } else
    {
      diff(x) < 0
    }
  )
  if(is_non_empty(bad_diff_x))
  {
    d <- data.frame(
      Position = paste(bad_diff_x, bad_diff_x + 1, sep = "/"),
      ValueBefore = x[bad_diff_x],  
      ValueAfter = x[bad_diff_x + 1]
    )
    msg <- gettext(
      if(strictly)
      {
        "The values of %s are not strictly monotonic increasing.\n%s"
      } else
      {
        "The values of %s are not monotonic increasing.\n%s"
      }
    )
    return(false(msg, .xname, assertive.base:::print_and_capture(d)))
  }
  TRUE
}

#' @rdname is_monotonic_increasing
#' @export
is_monotonic_decreasing <- function(x, strictly = FALSE, .xname = get_name_in_parent(x))
{
  bad_diff_x <- which(
    if(strictly)
    {
      diff(x) >= 0
    } else
    {
      diff(x) > 0
    }
  )
  if(is_non_empty(bad_diff_x))
  {
    d <- data.frame(
      Position = paste(bad_diff_x, bad_diff_x + 1, sep = "/"),
      ValueBefore = x[bad_diff_x],  
      ValueAfter = x[bad_diff_x + 1]
    )
    msg <- gettext(
      if(strictly)
      {
        "The values of %s are not strictly monotonic decreasing.\n%s"
      } else
      {
        "The values of %s are not monotonic decreasing.\n%s"
      }
    )
    return(false(msg, .xname, assertive.base:::print_and_capture(d)))
  }
  TRUE
}
