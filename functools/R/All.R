#' Test if all items in an object evaluate to TRUE.
#'
#' \code{All()} is a predicate functional that takes a predicate function
#' \code{.f} and an iterable object \code{.x} and:
#' \enumerate{
#'   \item iterates over each item \code{i} in object \code{.x},
#'   \item evaluates \code{.f(i)},
#'   \item and ultimately returns TRUE if all items \code{i} in object \code{.x} evaluate as TRUE.
#'   }
#'
#' @param .x An iterable object.
#' @param .f A predicate function.
#' @param ... Further arguments passed to the predicate function.
#' @param na.rm A logical value indicating whether NA values should be stripped before the computation proceeds.
#' @return A logical value indicating if all items evaluated as TRUE.
#' @family predicate functionals
#' @seealso \code{\link{Any}} to test if all items in an object evaluate to TRUE.
#' @examples
#' # Examples
#' data(mtcars)
#' All(mtcars, is.numeric) # TRUE
#' All(mtcars, is.character) # FALSE
#' mtcars$am <- factor(mtcars$am)
#' All(mtcars, is.numeric) # FALSE
#' All(mtcars, is.factor) # FALSE
#'
#' # Handles NAs and NULLs
#' All(list(NA, "3", NULL), is.numeric) # FALSE
#' All(list(NA, 3, NULL), is.numeric) # FALSE
#' All(list(NA, "3", NULL, 5), is.numeric) # FALSE
#'
#' # Use na.rm = TRUE to remove NAs and NULLS
#' All(list(NA, TRUE), Identity) # NA
#' All(list(NA, TRUE), Identity, na.rm = TRUE) # TRUE
#' @export
All <- function(.x, .f, ..., na.rm = FALSE) {
  .f <- match.fun(.f)
  return(all(vapply(.x, .f, logical(1), ...), na.rm = na.rm))
}
