#' Is the input a condition?
#'
#' Checks to see if the input is a message, warning or error.
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
#' # stop returns simpleErrors, unless wrapped in a call to try
#' simple_err <- tryCatch(stop("!!!"), error = function(e) e)
#' is_simple_error(simple_err)
#' try_err <- try(stop("!!!"))
#' is_try_error(try_err)
#' 
#' # is_error checks for both error types
#' is_error(try_err)
#' is_error(simple_err)
#' 
#' # warning returns simpleWarnings
#' simple_warn <- tryCatch(warning("!!!"), warning = function(w) w)
#' is_simple_warning(simple_warn)
#' is_warning(simple_warn)
#' 
#' # message returns simpleMessages
#' simple_msg <- tryCatch(message("!!!"), message = function(m) m)
#' is_simple_message(simple_msg)
#' is_message(simple_msg)
#' 
#' # These examples should fail.
#' assertive.base::dont_stop(assert_is_simple_error(try_err))
#' assertive.base::dont_stop(assert_is_try_error(simple_err))
#' @importFrom assertive.base is2
#' @export
is_try_error <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "try-error", .xname)
}

#' @rdname is_try_error
#' @export
is_simple_error <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "simpleError", .xname)
}

#' @rdname is_try_error
#' @export
is_error <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "error", .xname)
}

#' @rdname is_try_error
#' @export
is_simple_warning <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "simpleWarning", .xname)
}

#' @rdname is_try_error
#' @export
is_warning <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "warning", .xname)
}

#' @rdname is_try_error
#' @export
is_simple_message <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "simpleMessage", .xname)
}

#' @rdname is_try_error
#' @export
is_message <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "message", .xname)
}

#' @rdname is_try_error
#' @export
is_condition <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "condition", .xname)
}
