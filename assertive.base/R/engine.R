#' Throws an error if a condition isn't met
#'
#' The workhorse of the package.  If a condition isn't met, then an error
#' is thrown.  This function is exported for use by package developers so
#' that they can create their own assert functions.  
#'
#' @param predicate Function that returns a logical value (possibly 
#' a vector).
#' @param ... Passed to the \code{predicate} function.
#' @param msg The error message, in the event of failure.
#' @param what Either 'all' or 'any', to reduce vectorised tests to a 
#' single value.
#' @param na_ignore A logical value.  If \code{FALSE}, \code{NA} values
#' cause an error; otherwise they do not.  Like \code{na.rm} in many
#' stats package functions, except that the position of the failing
#' values does not change.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{FALSE} with the attribute \code{message}, as provided
#' in the input.
#' @note Missing values are considered as \code{FALSE} for the purposes of
#' whether or not an error is thrown.
#' @export
assert_engine <- function(predicate, ..., msg = "The assertion failed.", what = c("all", "any"), na_ignore = FALSE, severity = c("stop", "warning", "message", "none"))
{
  handler_type <- match.arg(severity)
  dots <- list(...)
  return_value <- if(length(dots) > 0) dots[[1]] else NULL
  if(handler_type == "none") 
  {
    return(invisible(return_value))
  }
  what <- match.fun(match.arg(what))
  predicate_name <- get_name_in_parent(predicate)
  
  ok <- predicate(...)
  if(inherits(ok, "scalar_with_cause"))
  {
    if(!isTRUE(ok))
    {
      if(missing(msg))
      {
        msg <- cause(ok)
      }
      give_feedback(handler_type, msg, predicate_name)
    }
  } else # inherits(ok, "vector_with_cause")
  {
    really_ok <- if(na_ignore)
    {
      # ok can be TRUE or NA; FALSE is bad
      ok | is.na(ok)
    } else
    {
      # ok can be TRUE; FALSE or NA is bad
      ok & !is.na(ok)
    }
    if(!what(really_ok))
    {
      # Append first few failure values and positions to the error message.
      msg <- paste(enc2utf8(msg), print_and_capture(ok), sep = "\n")
      give_feedback(handler_type, msg, predicate_name)
    }
  }
  invisible(return_value)
}

give_feedback <- function(handler_type, msg, predicate_name)
{
  handler <- match.fun(
    handler_type
  )
  ass_condition <- switch(
    handler_type,
    stop = assertionError,
    warning = assertionWarning,
    message = assertionMessage
  )
  # Throw error/warning/message
  caller <- sys.call(-3)
  # UTF-8 characters do not display correctly under Windows for some 
  # LC_CTYPE locale values, but there isn't much assertive can do about that.
  # http://stackoverflow.com/q/32696241/134830
  handler(ass_condition(paste(predicate_name, msg, sep = " : "), caller, predicate_name))
}

#' FALSE, with a cause of failure.
#'
#' Always returns the value \code{FALSE}, with a cause attribute.
#'
#' @param ... Passed to \code{gettextf} to create a cause of failure message.
#' @return \code{FALSE} with the attribute \code{cause}, as provided
#' in the input.
#' @seealso \code{\link{cause}} and \code{\link{na}}.
#' @export
false <- function(...)
{
  msg <- if(nargs() > 0L) sprintf(...) else ""
  x <- FALSE
  cause(x) <- msg[1]
  class(x) <- c("scalar_with_cause", "logical")
  x
}

#' NA, with a cause of failure.
#'
#' Always returns the value (logical) \code{NA}, with a cause attribute.
#'
#' @param ... Passed to \code{gettextf} to create a cause of failure message.
#' @return \code{NA} with the attribute \code{cause}, as provided
#' in the input.
#' @seealso \code{\link{cause}} and \code{\link{false}}.
#' @export
na <- function(...)
{
  msg <- if(nargs() > 0L) sprintf(...) else ""
  x <- NA
  cause(x) <- msg[1]
  class(x) <- c("scalar_with_cause", "logical")
  x
}
