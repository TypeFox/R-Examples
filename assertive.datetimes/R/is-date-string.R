#' Does the character vector contain dates? 
#' 
#' Checks that the input contains dates or times.
#' 
#' @param x Input to check.
#' @param format Expected format of the dates.  See 
#' \code{\link[base]{strptime}}.
#' @param na_ignore A logical value.  If \code{FALSE}, \code{NA} values
#' cause an error; otherwise they do not.  Like \code{na.rm} in many
#' stats package functions, except that the position of the failing
#' values does not change.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @param .xname Not intended to be used directly.
#' @return A logical vector that is \code{TRUE} when the input contains valid 
#' dates or times.
#' @examples
#' x <- c("9999-12-31 23:59:59", "wednesday", NA)
#' is_date_string(x)
#' assert_all_are_date_strings("01Aug1979", format = "%d%b%Y") #My DOB!
#' @seealso \code{\link[base]{strptime}} for specifying formats, and the 
#' \code{lubridate} package for automatic guessing of date formats (and other 
#' date manipulation functions).
#' @importFrom assertive.base use_first
#' @export
is_date_string <- function(x, format = "%F %T", .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "character", .xname)
  format <- use_first(format)
  call_and_name(
    function(x) 
    {
      ok <- !is.na(strptime(x, format))
      set_cause(ok, ifelse(is.na(x), "missing", "bad format"))
    }, 
    x
  )  
}

