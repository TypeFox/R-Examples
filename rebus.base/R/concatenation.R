#' Combine strings together
#'
#' Operator equivalent of \code{regex}.
#'
#' @param x A character vector.
#' @param y A character vector.
#' @return A character vector representing part or all of a regular expression.
#' @note \code{\%c\%} was the original operator for this ('c' for
#' 'concatenate').  This is hard work to type on a QWERTY keyboard
#' though, so it has been replaced with \code{\%R\%}.
#' @seealso \code{\link{regex}}, \code{\link[base]{paste}}
#' @examples
#' letters %R% LETTERS
#' @name Concatenation
#' @export
`%c%` <- function(x, y)
{
  .Deprecated("%R%")
  regex(x, y)
}

#' @rdname Concatenation
#' @export
`%R%` <- function(x, y)
{
  regex(x, y)
}
