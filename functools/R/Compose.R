#' Compose multiple functions.
#'
#' In infix and prefix forms.
#'
#' @param ... n functions to apply in order from right to left.
#' @param .f A function.
#' @param .g A function.
#' @return A function that will apply each function in order from right to left.
#' @family function operators
#' @examples
#' not_null <- `!` %O% is.null
#' not_null(4)
#' not_null(NULL)
#'
#' add1 <- function(x) x + 1
#' Compose(add1,add1)(8)
#' @export
Compose <- function(...) {
  fs <- lapply(list(...), match.fun)
  n <- length(fs)

  last <- fs[[n]]
  rest <- fs[-n]

  function(...) {
    out <- last(...)
    for (f in rev(rest)) {
      out <- f(out)
    }
    out
  }
}

#' @rdname Compose
#' @export
#' @usage .f \%O\% .g
'%O%' <- function(.f, .g) {
  .f <- match.fun(.f)
  .g <- match.fun(.g)
  function(...) {
    .f(.g(...))
  }
}
