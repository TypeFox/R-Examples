#' Character classes
#'
#' Match character classes.
#' @param lo A non-negative integer. Minimum number of repeats, when grouped.
#' @param hi positive integer. Maximum number of repeats, when grouped.
#' @param char_class A logical value. Should \code{x} be wrapped in a character
#' class?  If \code{NA}, the function guesses whether that's a good idea.
#' @return A character vector representing part or all of a regular expression.
#' @note R has many built-in locale-dependent character classes, like
#' \code{[:alnum:]} (representing lower or upper case letters or numbers).
#' There are also some generic character classes like \code{\\w} (representing
#' lower or upper case letters or numbers or underscores).
#' Finally, there are ASCII-only ways of specifying letters like \code{a-zA-Z}.
#' Which version you want depends upon how you want to deal with international
#' characters, and the vagaries of the underlying regular expression engine.
#' I suggest reading the \code{\link[base]{regex}} help page and doing lots of
#' testing.
#' @references \url{http://www.regular-expressions.info/shorthand.html} and
#' \url{http://www.rexegg.com/regex-quickstart.html#posix}
#' @seealso \code{\link[base]{regex}}, \code{\link[rebus.unicode]{Unicode}}
#' @examples
#' # R character classes
#' alnum()
#' alpha()
#' blank()
#' cntrl()
#' digit()
#' graph()
#' lower()
#' printable()
#' punct()
#' space()
#' upper()
#' hex_digit()
#'
#' # Generic classes
#' any_char()
#' dgt()
#' wrd()
#' spc()
#'
#' # Generic negated classes
#' not_dgt()
#' not_wrd()
#' not_spc()
#'
#' # Non-locale-specific classes
#' ascii_digit()
#' ascii_lower()
#' ascii_upper()
#'
#' # Don't provide a class wrapper
#' digit(char_class = FALSE) # same as DIGIT
#'
#' # Match repeated values
#' digit(3)
#' digit(3, 5)
#' digit(0)
#' digit(1)
#' digit(0, 1)
#'
#' # Ranges of characters
#' char_range(0, 7) # octal number
#'
#' # Usage
#' (rx <- digit(3))
#' stringi::stri_detect_regex(c("123", "one23"), rx)
#' @include constants.R
#' @include grouping-and-repetition.R
#' @name ClassGroups
NULL

#' @rdname ClassGroups
#' @export
alnum <- function(lo, hi, char_class = TRUE)
{
  repeated(ALNUM, lo, hi, char_class = char_class)
}

#' @rdname ClassGroups
#' @export
alpha <- function(lo, hi, char_class = TRUE)
{
  repeated(ALPHA, lo, hi, char_class = char_class)
}

#' @rdname ClassGroups
#' @export
blank <- function(lo, hi, char_class = TRUE)
{
  repeated(BLANK, lo, hi, char_class = char_class)
}

#' @rdname ClassGroups
#' @export
cntrl <- function(lo, hi, char_class = TRUE)
{
  repeated(CNTRL, lo, hi, char_class = char_class)
}

#' @rdname ClassGroups
#' @export
digit <- function(lo, hi, char_class = TRUE)
{
  repeated(DIGIT, lo, hi, char_class = char_class)
}

#' @rdname ClassGroups
#' @export
graph <- function(lo, hi, char_class = TRUE)
{
  repeated(GRAPH, lo, hi, char_class = char_class)
}

#' @rdname ClassGroups
#' @export
lower <- function(lo, hi, char_class = TRUE)
{
  repeated(LOWER, lo, hi, char_class = char_class)
}

#' @rdname ClassGroups
#' @export
printable <- function(lo, hi, char_class = TRUE)
{
  repeated(PRINT, lo, hi, char_class = char_class)
}

#' @rdname ClassGroups
#' @export
punct <- function(lo, hi, char_class = TRUE)
{
  repeated(PUNCT, lo, hi, char_class = char_class)
}

#' @rdname ClassGroups
#' @export
space <- function(lo, hi, char_class = TRUE)
{
  repeated(SPACE, lo, hi, char_class = char_class)
}

#' @rdname ClassGroups
#' @export
upper <- function(lo, hi, char_class = TRUE)
{
  repeated(UPPER, lo, hi, char_class = char_class)
}

#' @rdname ClassGroups
#' @export
hex_digit <- function(lo, hi, char_class = TRUE)
{
  repeated(HEX_DIGIT, lo, hi, char_class = char_class)
}

#' @rdname ClassGroups
#' @export
any_char <- function(lo, hi)
{
  repeated(ANY_CHAR, lo, hi, char_class = FALSE)
}

#' @rdname ClassGroups
#' @export
dgt <- function(lo, hi, char_class = TRUE)
{
  repeated(DGT, lo, hi, char_class = char_class)
}

#' @rdname ClassGroups
#' @export
wrd <- function(lo, hi, char_class = TRUE)
{
  repeated(WRD, lo, hi, char_class = char_class)
}

#' @rdname ClassGroups
#' @export
spc <- function(lo, hi, char_class = TRUE)
{
  repeated(SPC, lo, hi, char_class = char_class)
}

#' @rdname ClassGroups
#' @export
not_dgt <- function(lo, hi, char_class = TRUE)
{
  repeated(NOT_DGT, lo, hi, char_class = char_class)
}

#' @rdname ClassGroups
#' @export
not_wrd <- function(lo, hi, char_class = TRUE)
{
  repeated(NOT_WRD, lo, hi, char_class = char_class)
}

#' @rdname ClassGroups
#' @export
not_spc <- function(lo, hi, char_class = TRUE)
{
  repeated(NOT_SPC, lo, hi, char_class = char_class)
}

#' @rdname ClassGroups
#' @export
ascii_digit <- function(lo, hi, char_class = TRUE)
{
  repeated(ASCII_DIGIT, lo, hi, char_class = char_class)
}

#' @rdname ClassGroups
#' @export
ascii_lower <- function(lo, hi, char_class = TRUE)
{
  repeated(ASCII_LOWER, lo, hi, char_class = char_class)
}

#' @rdname ClassGroups
#' @export
ascii_upper <- function(lo, hi, char_class = TRUE)
{
  repeated(ASCII_UPPER, lo, hi, char_class = char_class)
}

#' @rdname ClassGroups
#' @export
ascii_alpha <- function(lo, hi, char_class = TRUE)
{
  repeated(ASCII_LOWER %R% ASCII_UPPER, lo, hi, char_class = char_class)
}

#' @rdname ClassGroups
#' @export
ascii_alnum <- function(lo, hi, char_class = TRUE)
{
  repeated(ASCII_ALPHA %R% ASCII_DIGIT, lo, hi, char_class = char_class)
}

#' @rdname ClassGroups
#' @export
char_range <- function(lo, hi, char_class = lo < hi)
{
  lo <- get_first_char(lo)
  hi <- get_first_char(hi)

  x <- if(lo < hi)
  {
    regex(lo, "-", hi)
  } else if(lo == hi)
  {
    warning("'lo' and 'hi' are the same value.  Return 'lo'.")
    as.regex(lo)
  } else # lo > hi
  {
    stop("'hi' is less than 'lo'.")
  }
  if(char_class)
  {
    x <- char_class(x)
  }
  x
}

get_first_char <- function(x)
{
  x <- as.character(x)
  if(any(is.na(x)))
  {
    stop("Missing values are not allowed.")
  }
  if(nchar(x) > 1)
  {
    warning("Returning only the first character from x.")
    return(substring(x, 1, 1))
  }
  x
}
