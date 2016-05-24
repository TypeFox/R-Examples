#' Create a function that that always returns a specific object.
#'
#' \code{Always(.x)} is a closure function that takes any object \code{.x},
#' and returns a function that always returns object \code{.x}.
#'
#' @param .x An object.
#' @return A function that itself returns \code{.x}.
#' @family closures
#' @examples
#' # comment here
#' always_0 <- Always(0)
#' always_0() # 0
#' always_true <- Always(TRUE)
#' always_true() # TRUE
#'
#' @export
Always <- function(.x) {
  return(function() {return(.x)})
}
