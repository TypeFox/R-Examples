#' Backreferences
#'
#' Backreferences for replacement operations.
#' @references \url{http://www.regular-expressions.info/backref.html} and
#' \url{http://www.rexegg.com/regex-capture.html}
#' @examples
#' # For R's PCRE and Perl engines
#' REF1
#' REF2
#' # and so on, up to
#' REF9
#'
#' # For stringi/stringr's ICU engine
#' ICU_REF1
#' ICU_REF2
#' # and so on, up to
#' ICU_REF9
#'
#' # Usage
#' sub("a(b)c(d)", REF1 %R% REF2, "abcd")
#' stringi::stri_replace_first_regex("abcd", "a(b)c(d)", ICU_REF1 %R% ICU_REF2)
#' @name Backreferences
#' @include regex-methods.R
#' @export
REF1 <- regex("\\1")

#' @name Backreferences
#' @export
REF2 <- regex("\\2")

#' @name Backreferences
#' @export
REF3 <- regex("\\3")

#' @name Backreferences
#' @export
REF4 <- regex("\\4")

#' @name Backreferences
#' @export
REF5 <- regex("\\5")

#' @name Backreferences
#' @export
REF6 <- regex("\\6")

#' @name Backreferences
#' @export
REF7 <- regex("\\7")

#' @name Backreferences
#' @export
REF8 <- regex("\\8")

#' @name Backreferences
#' @export
REF9 <- regex("\\9")

#' @name Backreferences
#' @export
ICU_REF1 <- regex("$1")

#' @name Backreferences
#' @export
ICU_REF2 <- regex("$2")

#' @name Backreferences
#' @export
ICU_REF3 <- regex("$3")

#' @name Backreferences
#' @export
ICU_REF4 <- regex("$4")

#' @name Backreferences
#' @export
ICU_REF5 <- regex("$5")

#' @name Backreferences
#' @export
ICU_REF6 <- regex("$6")

#' @name Backreferences
#' @export
ICU_REF7 <- regex("$7")

#' @name Backreferences
#' @export
ICU_REF8 <- regex("$8")

#' @name Backreferences
#' @export
ICU_REF9 <- regex("$9")


#' Make the regular expression recursive.
#'
#' Makes the regular expression (or part of it) recursive.
#' @param x A character vector.
#' @return A character vector representing part or all of a regular expression.
#' @note Recursion is not supported by R's internal PRCE engine or
#' \code{stringi}'s ICU engine.
#' @references \url{http://www.regular-expressions.info/recurse.html} and
#' \url{http://www.rexegg.com/regex-recursion.html}
#' @examples
#' recursive("a")
#'
#' # Recursion isn't supported by R's PRCE engine or
#' # stringi/stringr's ICU engine
#' x <- c("ab222z", "ababz", "ab", "abab")
#' rx <- "ab(?R)?z"
#' grepl(rx, x, perl = TRUE)
#' try(grepl(rx, x))
#' try(stringi::stri_detect_regex(x, rx))
#' @export
recursive <- function(x)
{
  regex(x, "(?R)")
}
