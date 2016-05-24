#' Is the input TRUE/FALSE/NA?
#' 
#' Checks to see if the input is \code{TRUE},  \code{FALSE} or  \code{NA}.
#'
#' @param x Input to check.  See note.
#' @param allow_attributes If \code{TRUE}, a scalar value of \code{TRUE}
#' with attributes is allowed.
#' @param coerce_to_logical Logical: should the input be coerced to logical
#' before checking?  See note.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @note \code{is_identical_to_true} wraps the base function \code{isTRUE}, 
#' providing more information on failure.  Likewise, 
#' \code{is_identical_to_false} checks that the input is identical to FALSE.  If
#' \code{allow_attributes} is \code{TRUE}, a scalar value of \code{TRUE} with 
#' attributes is allowed. \code{is_true} and \code{is_false} are vectorized, 
#' returning \code{TRUE} when the inputs are \code{TRUE} and \code{FALSE} 
#' respectively.
#' 
#' The for \code{is_true}, \code{is_false}, \code{is_not_true} and 
#' \code{is_not_false}, \code{x} argument will be coerced to be a logical vector 
#' if it isn't already.  
#' 
#' Coercion to logical is optional for \code{is_na} and \code{is_not_na}. If 
#' you do coerce, it means that \code{is_na} differs in behaviour from 
#' \code{base::is.na} for character vector, list and data frame inputs.  To 
#' replicate the behaviour of \code{is.na}, ensure the argument 
#' \code{coerce_to_logical} is \code{FALSE} (this is the default).
#' 
#' Note that in assertive version 0.1-4 and prior, 
#' \code{is_identical_to_true/false} were named \code{is_true/false} and the 
#' vectorized versions were not present.
#' @return The \code{is*} functions return \code{TRUE} if the input is 
#' \code{TRUE}/\code{FALSE}. The \code{assert_*} functions return nothing but 
#' throw an error if the corresponding \code{is_*} function returns 
#' \code{FALSE}.
#' @seealso \code{\link[base]{isTRUE}}.
#' @examples
#' # Checks against logical values using base::identical
#' assert_is_identical_to_true(TRUE)
#' assert_is_identical_to_false(FALSE)
#' assert_is_identical_to_na(NA)
#' 
#' # Other NA types match
#' assert_is_identical_to_na(NA_complex_)
#' 
#' # NaN is not NA
#' dont_stop(assert_is_identical_to_na(NaN))
#' 
#' # For a slightly less strict test, you can ignore attributes
#' assert_is_identical_to_true(c(truth = TRUE), allow_attributes = TRUE)
#' assert_is_identical_to_false(matrix(FALSE), allow_attributes = TRUE)
#' assert_is_identical_to_na(structure(NA, class = "nanana"), allow_attributes = TRUE)
#' 
#' # Vectorized predicates (package name explicitly given to prevent
#' # problems with testthat name clash)
#' x <- c(TRUE, FALSE, NA)
#' assertive.base::is_true(x)
#' assertive.base::is_false(x)
#' is_na(x)
#' 
#' # ...and their opposites
#' is_not_true(x)
#' is_not_false(x)
#' is_not_na(x)
#' 
#' # Check that at least one element fits the condition
#' assert_any_are_true(x)
#' assert_any_are_false(x)
#' assert_any_are_na(x)
#' 
#' # These checks should fail:
#' dont_stop({
#'   assert_is_identical_to_true(c(truth = TRUE))
#'   assert_is_identical_to_true(1)
#'   assert_is_identical_to_true(c(TRUE, TRUE))
#'   assert_is_identical_to_false(matrix(FALSE))
#'   assert_is_identical_to_na(structure(NA, class = "nanana"))
#'   assert_all_are_true(x)
#'   assert_all_are_false(x)
#'   assert_all_are_na(x)
#' })
#' 
#' # base::is.na has non-standard behaviour for data.frames and lists.
#' # is_na and is_not_na coerce to logical vectors (except character input).
#' # unlist the input or use an apply function.
#' d <- data.frame(
#'   x = c(TRUE, FALSE, NA), 
#'   y = c(0, NA, 2), 
#'   z = c("a", "NA", NA)
#' )
#' is.na(d)
#' is_na(unlist(d))
#' @name Truth
NULL

#' @rdname Truth
#' @export
assert_is_identical_to_false <- function(x, allow_attributes = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                  
  assert_engine(
    is_identical_to_false,
    x, 
    allow_attributes = allow_attributes, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname Truth
#' @export
assert_is_identical_to_na <- function(x, allow_attributes = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                  
  assert_engine(
    is_identical_to_na,
    x,
    allow_attributes = allow_attributes, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname Truth
#' @export
assert_is_identical_to_true <- function(x, allow_attributes = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                  
  assert_engine(
    is_identical_to_true,
    x,
    allow_attributes = allow_attributes, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}
