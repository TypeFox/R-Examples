#' ISO 8601 date-time classes
#'
#' Match ISO 8601 date and time classes.
#' @param lo A non-negative integer. Minimum number of repeats, when grouped.
#' @param hi positive integer. Maximum number of repeats, when grouped.
#' @param char_class \code{TRUE} or \code{FALSE}. Should the values be wrapped
#' into a character class?
#' @return A character vector representing part or all of a regular expression.
#' @references \url{http://www.iso.org/iso/iso8601}
#' @examples
#' iso_date()
#' iso_time()
#' iso_datetime()
#'
#' # With repetition
#' iso_date(3, 6)
#' iso_time(1, Inf)
#' iso_datetime(0, Inf)
#'
#' # Without a class wrapper
#' iso_date(char_class = FALSE)
#' @name IsoClasses
#' @aliases IsoDateTime
#' @export
iso_date <- function(lo, hi, char_class = TRUE)
{
  repeated(ISO_DATE, lo, hi, char_class = char_class)
}

#' @rdname IsoClasses
#' @export
iso_time <- function(lo, hi, char_class = TRUE)
{
  repeated(ISO_TIME, lo, hi, char_class = char_class)
}

#' @rdname IsoClasses
#' @export
iso_datetime <- function(lo, hi, char_class = TRUE)
{
  repeated(ISO_DATETIME, lo, hi, char_class = char_class)
}
