#' Get or set the \code{"cause"} attribute
#'
#' Gets or sets the \code{"cause"} (of failure) attribute of a variable.
#'
#' @param x Any variable.
#' @param value Passed to \code{gettextf} and stored in the \code{"cause"}
#' attribute.
#' @return The get method returns the \code{"cause"} attribute.
#' @seealso \code{\link{set_cause}}
#' @examples
#' # Scalar case
#' yn <- is_identical_to_true(FALSE)
#' cause(yn)
#' 
#' # Vector case
#' yn <- is_true(c(TRUE, FALSE, NA))
#' cause(yn)
#' @export
cause <- function(x)
{
  y <- attr(x, "cause")
  if(is.null(y))
  {
    return(noquote(character(length(x))))
  }
  y
}

#' @rdname cause
#' @export
`cause<-` <- function(x, value)
{
  # Can't use is_scalar here due to dependency on this
  if(length(value) != 1 && length(value) != length(x))
  {
    stop(
      sprintf(
        "The length of value should be 1 or the length of x (%d) but is %d.", 
        length(x),
        length(value)
      )
    )
  }
  attr(x, "cause") <- noquote(as.character(value))
  x
}

#' Set a cause and return the input
#' 
#' Sets the cause attribute of an object and returns that object.
#' @param x A variable.
#' @param false_value A character vector to set the cause to, where \code{x} is
#' \code{FALSE}.
#' @param missing_value A character vector to set the cause to, where \code{x} is
#' \code{NA}.
#' @details If \code{x} is \code{TRUE} everywhere, this returns the input 
#' without setting a cause.  Otherwise, the cause is an empty string where 
#' \code{x} is \code{TRUE}, \code{false_value} where it is \code{FALSE}, and
#' \code{missing_value} where it is \code{NA}.
#' @return \code{x}, with a new cause attribute.
#' @seealso \code{\link{cause}}, \code{\link[stats]{setNames}}
#' @export
set_cause <- function(x, false_value, missing_value = "missing")
{
  if(all(!is.na(x) & x)) return(x)
  len_x <- length(x)
  # TRUEs
  cause_value <- character(len_x)
  # NAs
  missing_value <- rep_len(missing_value, len_x)
  missing_index <- is.na(x)
  cause_value[missing_index] <- missing_value[missing_index]
  # FALSEs
  false_value <- rep_len(false_value, len_x)
  false_index <- !x & !is.na(x)
  cause_value[false_index] <- false_value[false_index]
  
  cause(x) <- cause_value
  class(x) <- c("vector_with_cause", "logical")
  x
}

#' @rdname print.vector_with_cause
#' @method print scalar_with_cause
#' @export
print.scalar_with_cause <- function(x, ...)
{
  if(length(x) != 1L)
  {
    stop("x is malformed; it should have length 1.", domain = NA) 
  }
  print(x[1])
  cat("Cause of failure: ", cause(x), "\n")
}    

#' Print methods for objects with a cause attribute
#' 
#' Prints objects of class \code{scalar_with_cause} and 
#' \code{vector_with_cause}.
#' @param x an object of class \code{scalar_with_cause} or
#' \code{vector_with_cause}.
#' @param na_ignore A logical value.  If \code{FALSE}, \code{NA} values
#' are printed; otherwise they do not.  Like \code{na.rm} in many
#' stats package functions, except that the position of the failing
#' values does not change.
#' @param n_to_show A natural number.  The maximum number of failures 
#' to show. 
#' @param ... Currently unused.
#' @method print vector_with_cause
#' @importFrom utils head
#' @export
print.vector_with_cause <- function(x, na_ignore = FALSE, n_to_show = 10, ...)
{
  cause_x <- cause(x)
  names_x <- names(x)
  if(is.null(names_x))
  {
    names_x <- character(length(x))
  }
  x <- strip_attributes(x)
  ok <- if(na_ignore)
  {
    # ok can be TRUE or NA; FALSE is bad
    x | is.na(x)
  } else
  {
    # ok can be TRUE; FALSE or NA is bad
    x & !is.na(x)
  }
  
  # Append first few failure values and positions to the error message.
  fail_index <- which(!ok)
  n <- length(fail_index)
  fail_index <- head(fail_index, n_to_show)
  failures <- data.frame(
    Position = fail_index,
    Value    = truncate(names_x[fail_index]),
    Cause    = unclass(cause_x[fail_index]), # See bug 15997
    row.names = seq_along(fail_index)
  )
  # Slightly convoluted way of creating message done to ensure that xgettext
  # creates all the translation strings
  msg_showing_first <- if(nrow(failures) < n) 
  {
    paste0(
      " ", 
      gettextf(
        "(showing the first %d)", 
        nrow(failures), 
        domain = "R-assertive.base"
      )
    )
  } else ""
  msg_n_failures <- ngettext(
    n,
    "There was %d failure%s:\n",
    "There were %d failures%s:\n",
    domain = "R-assertive.base"
  )
  cat(enc2utf8(sprintf(msg_n_failures, n, msg_showing_first)))
  print(failures)
}
