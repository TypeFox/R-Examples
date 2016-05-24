first <- function(...) head(..., n = 1)
but_last <- function(...) head(..., n = -1)
last <- function(...) tail(..., n = 1)

length_zero <- function(x) length(x) == 0
length_one <- function(x) length(x) == 1

is_list_of <- function(xs, classes) {
  assert_that(is.list(xs))
  all(vapply(xs, function(x) inherits(x, classes), logical(1)))
}

merge_lists <- function(x, y) {
  x[names(y)] <- y
  x
}

#' Higher-order functions for dealing with lists
#'
#' These functions were inspired by underscore.js.
#'
#' @name list_functions
#' @param key the name of a value in a list
#' @param keys a character vector of names in a list
#' @param xss a list of lists
#' @return \code{pluck} returns an unnamed value and \code{pluck_apply} returns
#'   a list of unnamed values. \code{pick} returns a simplified version of the
#'   original list. \code{pick_apply} returns a list of simplified lists.
#'
#' @details \itemize{ \item \code{pluck}: Pluck a named value from a list \item
#'   \code{pick}: Simplify a list by picking out whitelisted names}
#'
#'   The simple versions of \code{pluck} and \code{pick} are curried functions,
#'   meaning that they return a function which can be applied to a list. See the
#'   syntax in the usage section.
#' @keywords internal
NULL

#' @rdname list_functions
pluck <- function(key) {
  function(xs) xs[[key]]
}

#' @rdname list_functions
pluck_apply <- function(key, xss) {
  assert_that(is_list_of(xss, "list"))
  lapply(xss, pluck(key))
}

#' @rdname list_functions
pick <- function(keys) {
  function(xs) {
    classes <- class(xs)
    xs <- xs[is.element(names(xs), keys)]
    class(xs) <- classes
    xs
  }
}

#' @rdname list_functions
pick_apply <- function(keys, xss) {
  assert_that(is_list_of(xss, "list"))
  classes <- class(xss)
  xss <- lapply(xss, pick(keys))
  class(xss) <- classes
  xss
}


