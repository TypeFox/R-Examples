#' @include workaround.R
#' @rdname is_in_past   
#' @export
is_after <- function(x, y, .xname = get_name_in_parent(x), 
  .yname = get_name_in_parent(y))
{
  if(!is_date(x)) 
  {
    x <- coerce_to(x, "POSIXct", .xname)
  }
  y <- if(is_date(x))
  {
    coerce_to(y, "Date", .yname)  
  } else
  {
    coerce_to(y, "POSIXct", .yname)
  }
  call_and_name_retro(
    function(x)
    {
      ok <- x > y
      set_cause(ok, "before")
    },
    x
  )
}

#' @rdname is_in_past   
#' @export
is_before <- function(x, y, .xname = get_name_in_parent(x), 
  .yname = get_name_in_parent(y))
{  
  if(!is_date(x)) 
  {
    x <- coerce_to(x, "POSIXct", .xname)
  }
  y <- if(is_date(x))
  {
    coerce_to(y, "Date", .yname)  
  } else
  {
    coerce_to(y, "POSIXct", .yname)
  }
  call_and_name_retro(
    function(x)
    {
      ok <- x < y
      set_cause(ok, "after")
    },
    x
  )
}

#' @rdname is_in_past   
#' @export
is_in_future <- function(x, .xname = get_name_in_parent(x))
{
  # Dates stay as dates, otherwise convert to POSIXct
  if(!is_date(x)) 
  {
    x <- coerce_to(x, "POSIXct", .xname)
  }
  call_and_name_retro(
    function(x)
    {
      ok <- x > Sys.time()
      set_cause(ok, "in past")
    },
    x
  )
}

#' Is the input in the past/future?
#'
#' Checks to see if the input is a time in the past/future, or before/after
#' some time point.
#'
#' @param x \code{Date} or \code{POSIXt} input to check.
#' @param y Another date-time object to compare against.
#' @param na_ignore A logical value.  If \code{FALSE}, \code{NA} values
#' cause an error; otherwise they do not.  Like \code{na.rm} in many
#' stats package functions, except that the position of the failing
#' values does not change.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @param .xname Not intended to be used directly.
#' @param .yname Not intended to be used directly.
#' @return The \code{is_*} function return \code{TRUE} if the input is 
#' a time in the future/past.  The \code{assert_*} functions return nothing but
#' throw an error if the corresponding \code{is_*} function returns 
#' \code{FALSE}.
#' @details The current time is determined by \code{Sys.time}, and the 
#' input is coerced to \code{POSIXct} format if necessary.
#' @note Note that the print method for \code{POSIXct} objects means that the
#' cause attribute (in the event of failures) is not shown.  You can still 
#' access it via, e.g., \code{cause(is_in_past(x))}.
#' @seealso \code{\link{Sys.time}}.
#' @examples
#' x <- Sys.time() + c(-1, 100)
#' is_in_past(x)
#' is_in_future(x)
#' 
#' # more generally, compare against any date-time
#' is_before(x, as.POSIXct("9999-12-31"))
#' is_after(x, as.POSIXct("0001-01-01"))
#' @importFrom assertive.types is_date
#' @export
is_in_past <- function(x, .xname = get_name_in_parent(x))
{
  # Dates stay as dates, otherwise convert to POSIXct
  if(!is_date(x)) 
  {
    x <- coerce_to(x, "POSIXct", .xname)
  }
  call_and_name_retro(
    function(x)
    {
      ok <- x < Sys.time()
      set_cause(ok, "in future")
    },
    x
  )
}
