#' Alternation
#'
#' Match one string or another.
#'
#' @param x A character vector.
#' @param y A character vector.
#' @param ... Character vectors.
#' @param capture A logical value indicating whether or not the result should be
#' captured.  See note.
#' @return A character vector representing part or all of a regular expression.
#' @note \code{or} takes multiple character vector inputs and returns a
#' character vector of the inputs separated by pipes. \code{\%|\%} is an operator
#' interface to this function. \code{or1} takes a single character vector and
#' returns a string collapsed by pipes.
#'
#' When \code{capture} is \code{TRUE}, the values are wrapped in a capture
#' group (see \code{\link{capture}}).  When \code{capture} is \code{FALSE} (the
#' default for \code{or} and \code{or1}), the values are wrapped in a
#' non-capture group (see \code{\link{token}}).  When \code{capture} is
#' \code{NA}, (the case for \code{\%|\%}) the values are not wrapped in
#' anything.
#' @seealso \code{\link[base]{paste}}
#' @references \url{http://www.regular-expressions.info/alternation.html}
#' @examples
#' # or takes an arbitrary number of arguments and groups them without capture
#' or(letters, LETTERS, "foo")
#'
#' # or1 takes a single character vector
#' or1(c(letters, LETTERS, "foo")) # Not the same as before!
#'
#' # Capture the group
#' or1(letters, capture = TRUE)
#'
#' # Don't create a group
#' or1(letters, capture = NA)
#'
#' # The pipe operator doesn't group
#' letters %|% LETTERS %|% "foo"
#'
#' # Usage
#' (rx <- or("dog", "cat", "hippopotamus"))
#' stringi::stri_detect_regex(c("boondoggle", "caterwaul", "water-horse"), rx)
#' @export
or <- function(..., capture = FALSE)
{
  n_dots <- length(list(...))
  if(n_dots < 2)
  {
    warning(
      "'or' is intended to be called with at least 2 arguments in '...'. ",
      sprintf(ngettext(n_dots, "%d was passed.", "%d were passed."), n_dots),
      " Maybe you wanted 'or1' instead?"
    )
  }
  engroup(paste(..., sep = "|", collapse = NULL), capture)
}

#' @rdname or
#' @export
`%|%` <- function(x, y)
{
  or(x, y, capture = NA_character_)
}

#' @rdname or
#' @export
or1 <- function(x, capture = FALSE)
{
  engroup(paste0(x, collapse = "|"), capture)
}
