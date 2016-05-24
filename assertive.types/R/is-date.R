#' Is the input a date?
#'
#' Checks to see if the input is a Date or POSIXt object.
#'
#' @param x Input to check.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @param .xname Not intended to be used directly.
#' @return The \code{is_*} functions return \code{TRUE} or \code{FALSE} 
#' depending upon whether or not the input is a datetime object.
#' 
#' The \code{assert_*} functions return nothing but throw an error if the 
#' corresponding \code{is_*} function returns \code{FALSE}.
#' @examples
#' is_date(Sys.Date())
#' is_posixct(Sys.time())
#' 
#' # These examples should fail.
#' assertive.base::dont_stop(assert_is_date(Sys.time()))
#' @importFrom assertive.base is2
#' @export
is_date <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "Date", .xname)
}

#' @rdname is_date
#' @export
is_posixct <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "POSIXct", .xname)
}

#' @rdname is_date
#' @export
is_posixlt <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "POSIXlt", .xname)
}
