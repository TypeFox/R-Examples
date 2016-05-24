#' Fail with a default value.
#'
#' \code{Fail_With()} turns a function that throws an error into a function
#' that returns a default value when there is an error. The essence of
#' \code{Fail_With()} is simple: it is just a wrapper around \code{try()},
#' the function that captures errors and allows execution to continue.
#'
#' @param .default default value.
#' @param .f any function that throws an error.
#' @param .silent logical: should the report of error messages be suppressed?
#' @return a function that returns a default value when there's an error.
#' @family function operators
#' @export
Fail_With <- function(.default = NULL, .f, .silent = FALSE) {
  .f <- match.fun(.f)
  return(function(...) {
    out <- .default
    try(out <- .f(...), silent = .silent)
    return(out)
  })
}
