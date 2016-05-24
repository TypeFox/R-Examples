#' Truthy
#'
#' \code{Truthy()} returns TRUE or FALSE if an object is TRUE or not.
#' An object is is "TRUE" if it is exists and is TRUE.
#'
#' @param .x an object.
#' @return a logical value.
#' @family predicate functions
#' @examples
#' # Returns if a value exists or not:
#' Truthy(TRUE) # TRUE
#' Truthy(FALSE) # FALSE
#' Truthy(NULL) # FALSE
#' Truthy(NA) # FALSE
#' Truthy(2L) # TRUE
#' Truthy(1L) # TRUE
#' Truthy(0L) # FALSE
#' Truthy("a") # TRUE
#'
#' @export
Truthy <- function(.x) {
  return(Existy(.x) && identical(.x, True()))
}
