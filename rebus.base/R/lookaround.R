#' Lookaround
#'
#' Zero length matching.
#' @param x A character vector.
#' @return A character vector representing part or all of a regular expression.
#' @references \url{http://www.regular-expressions.info/lookaround.html} and
#' \url{http://www.rexegg.com/regex-lookarounds.html}
#' @note Lookbehind is not supported by R's PRCE engine.  Use R's Perl engine
#' or \code{stringi}/\code{stringr}'s ICU engine.
#' @examples
#' x <- "foo"
#' lookahead(x)
#' negative_lookahead(x)
#' lookbehind(x)
#' negative_lookbehind(x)
#'
#' # Usage
#' x <- c("mozambique", "qatar", "iraq")
#' (rx_neg_class <- "q" %R% negated_char_class("u"))
#' (rx_neg_lookahead <- "q" %R% negative_lookahead("u"))
#' stringi::stri_detect_regex(x, rx_neg_class)
#' stringi::stri_detect_regex(x, rx_neg_lookahead)
#'
#' # PRCE engine doesn't support lookbehind
#' x2 <- c("queen", "vacuum")
#' (rx_lookbehind <- lookbehind("q")) %R% "u"
#' stringi::stri_detect_regex(x2, rx_lookbehind)
#' try(grepl(rx_lookbehind, x2))
#' try(grepl(rx_lookbehind, x2, perl = TRUE))
#' @export
lookahead <- function(x)
{
  regex("(?=", x, ")")
}

#' @rdname lookahead
#' @export
negative_lookahead <- function(x)
{
  regex("(?!", x, ")")
}

#' @rdname lookahead
#' @export
lookbehind <- function(x)
{
  regex("(?<=", x, ")")
}

#' @rdname lookahead
#' @export
negative_lookbehind <- function(x)
{
  regex("(?<!", x, ")")
}
