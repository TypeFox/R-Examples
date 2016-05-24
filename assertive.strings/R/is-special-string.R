# TODO: rename as is_numeric_character for consistency?

#' Does the string contain a number/logical value?
#' 
#' Check to see if a character vector contains numeric/logical strings.
#' @param x A character vector.
#' @param na_ignore A logical value.  If \code{FALSE}, \code{NA} values
#' cause an error; otherwise they do not.  Like \code{na.rm} in many
#' stats package functions, except that the position of the failing
#' values does not change.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @param .xname Not intended to be used directly.
#' @return \code{is_numeric_string} returns a logical vector that is \code{TRUE}
#' when the string contains numbers.  The corresponding \code{assert_*} 
#' functions return nothing but throw an error on failure.
#' @examples
#' is_numeric_string(c("1", "1.1", "-1.1e1", "one", NA))
#' # R only treats certain capitalizations of "true" and "false" as logical
#' x <- c(
#'   "TRUE", "FALSE", "true", "false", "True", "False", "trUE", "FaLsE", 
#'   "T", "F", "t", "f"
#' )
#' is_logical_string(x)
#' 
#' assert_all_are_numeric_strings(c("1", "2.3", "-4.5", "6e7", "8E-9"))
#' assert_any_are_numeric_strings(c("1", "Not a number"))
#' @importFrom assertive.base is_not_na
#' @export
is_numeric_string <- function(x, .xname)
{
  x <- coerce_to(x, "character", .xname)
  ok <- call_and_name(
    function(x)
    {
      suppressWarnings(
        {
          numx <- as.numeric(x)
          is_not_na(numx)
        }
      )
    },
    x
  )
  set_cause(ok, ifelse(is.na(x), "missing", "bad format"))
}

#' @rdname is_numeric_string
#' @export
is_logical_string <- function(x, .xname)
{
  x <- coerce_to(x, "character", .xname)
  ok <- call_and_name(
    function(x)
    {
      suppressWarnings(
        {
          logx <- as.logical(x)
          is_not_na(logx)
        }
      )
    },
    x
  )
  set_cause(ok, ifelse(is.na(x), "missing", "bad format"))
}

#' Is the input a single character?
#' 
#' Checks to see if he unput is a single character.
#' @param x A character vector.
#' @param na_ignore A logical value.  If \code{FALSE}, \code{NA} values
#' cause an error; otherwise they do not.  Like \code{na.rm} in many
#' stats package functions, except that the position of the failing
#' values does not change.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @param .xname Not intended to be used directly.
#' @return \code{is_single_character} returns \code{TRUE} when the input is a 
#' single character (as determined by \code{nchar}; this excludes \code{NA}).
#' The \code{assert_*} functions return nothing but throw an error if the 
#' corresponding \code{is_*} function returns \code{FALSE}.
#' @note The behaviour of this function has changed so that \code{NA} inputs 
#' now return \code{NA} rather than \code{FALSE}.
#' @seealso \code{\link[base]{nchar}}
#' @examples
#' x <- c("", "a", "aa", NA)
#' is_single_character(x)
#' @export
is_single_character <- function(x, .xname)
{
  x <- coerce_to(x, "character", .xname)
  ok <- call_and_name(
    function(x)
    {
      nch_is_1 <- nchar(x) == 1
      is_na_x <- is.na(x)
      if(any(is_na_x))
      {
        message("New behaviour: NA inputs now return NA.")
        nch_is_1[is_na_x] <- NA
      }
      nch_is_1
    },
    x
  )
  set_cause(ok, ifelse(is.na(x), "missing", "bad format"))  
}
