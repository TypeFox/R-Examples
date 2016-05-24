#' Condition classes
#' 
#' Error, warning, and message classes derived from their simple equivalents.
#' @param message A string describing the problem.
#' @param call A call describing the source of the condition.
#' @param predicate_name A string naming the predicate that was called when the 
#' condition occured.
#' @return An object of class \code{assertionError}, \code{assertionWarning}, or
#' \code{assertionMessage}.
#' @note These objects behave the same as the standard-issue \code{simpleError},
#' \code{simpleWarning}, and \code{simpleMessage} objects from base-R.  The
#' extra class allows you to provide custom handling for assertions inside 
#' \code{tryCatch}.
#' @examples 
#' tryCatch(
#'   assert_all_are_true(FALSE), 
#'   error = function(e) 
#'   {
#'     if(inherits(e, "assertionCondition"))
#'     {
#'       # Handle assertions
#'       message("This is an assertion condition.")
#'       
#'       # Handle assertions cause by a specific predicate
#'       if(e$predicate_name == "is_true")
#'       {
#'       }
#'     } else
#'     {
#'       # Handle other error types
#'     }
#'   }
#' )
#' @export
assertionError <- function(message, call = NULL, predicate_name = NULL)
{
  class <- c("assertionError", "assertionCondition", "simpleError", "error", "condition")
  structure(
    list(
      message = as.character(message), 
      call = call,
      predicate_name = predicate_name
    ), 
    class = class
  )
}

#' @rdname assertionError
#' @export
assertionWarning <- function(message, call = NULL, predicate_name = NULL)
{
  class <- c("assertionWarning", "assertionCondition", "simpleWarning", "warning", "condition")
  structure(
    list(
      message = as.character(message), 
      call = call,
      predicate_name = predicate_name
    ), 
    class = class
  )
}

#' @rdname assertionError
#' @export
assertionMessage <- function(message, call = NULL, predicate_name = NULL)
{
  class <- c("assertionMessage", "assertionCondition", "simpleMessage", "message", "condition")
  structure(
    list(
      message = as.character(message), 
      call = call,
      predicate_name = predicate_name
    ), 
    class = class
  )
}

