#' Binary Infix Concatenation
#'
#' \describe{
#' \item{\code{\%+\%}:}{binary infix operator for strings.}
#' \item{\code{\% + \%}:}{like \code{\%+\%} but with a space between its
#' inputs.}
#' \item{\code{\%,\%}, \code{\%or\%}, \code{\%and\%}:}{infix versions of
#'       \code{\link{cc}}, \code{\link{cc_or}}, \code{\link{cc_and}}.}
#' }
#' @param x,y Character vectors.
#' @examples
#' v <- "important value"
#' v %+% "!"
#'
#' message("Two" % + % "words")
#' @name infix
NULL

#' @rdname infix
#' @export
`%+%` <- function(x, y) {
  paste0(x, y)
}

#' @rdname infix
#' @export
`% + %` <- function(x, y) {
  paste(x, y)
}

#' @rdname infix
#' @export
`%,%` <- function(x, y) {
  cc(c(x, y))
}

#' @rdname infix
#' @export
`%or%` <- function(x, y) {
  cc_or(c(x, y), oxford = FALSE)
}

#' @rdname infix
#' @export
`%and%` <- function(x, y) {
  cc_and(c(x, y), oxford = FALSE)
}
