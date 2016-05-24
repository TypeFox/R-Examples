#' Capture a token, or not
#'
#' Create a token to capture or not.
#' @param x A character vector.
#' @return A character vector representing part or all of a regular expression.
#' @references \url{http://www.regular-expressions.info/brackets.html}
#' @seealso \code{\link{or}} for more examples
#' @examples
#' x <- "foo"
#' capture(x)
#' group(x)
#'
#' # Usage
#' # capture is good with match functions
#' (rx_price <- capture(DOLLAR %R% digit(1, Inf) %R% DOT %R% digit(2)))
#' (rx_quantity <- capture(digit(1, Inf)))
#' (rx_all <- rx_price %R% " for " %R% rx_quantity)
#' stringi::stri_match_first_regex("The price was $123.99 for 12.", rx_all)
#'
#' # group is mostly used with alternation.  See ?or.
#' @export
capture <- function(x)
{
  regex("(", x, ")")
}

#' @rdname capture
#' @export
group <- function(x)
{
  regex("(?:", x, ")")
}

#' @rdname capture
#' @export
token <- function(x)
{
  .Deprecated("group")
  group(x)
}

#' @rdname capture
#' @param capture Logical If \code{TRUE}, call \code{capture}; if \code{FALSE},
#' call \code{group}.
#' @export
engroup <- function(x, capture)
{
  switch(
    as.character(capture),
    "TRUE"  = capture(x),
    "FALSE" = group(x),
    as.regex(x)
  )
}
