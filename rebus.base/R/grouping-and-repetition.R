#' A range or char_class of characters.
#'
#' Match a range or char_class of characters.
#' @param ... Character vectors.
#' @return A character vector representing part or all of a regular expression.
#' @references \url{http://www.regular-expressions.info/charclass.html}
#' @examples
#' char_class(LOWER, "._")
#' negated_char_class(LOWER, "._")
#'
#' # Usage
#' (rx_odd <- char_class(1, 3, 5, 7, 9))
#' stringi::stri_detect_regex((1:10) ^ 2, rx_odd)
#' @export
char_class <- function(...)
{
  regex("[", ..., "]")
}

#' @rdname char_class
#' @export
negated_char_class <- function(...)
{
  regex("[^", ..., "]")
}

#' @rdname char_class
#' @export
negate_and_group <- function(...)
{
  .Defunct("negated_char_class")
  negated_char_class(...)
}

#' Repeat values
#'
#' Match repeated values.
#' @param x A character vector.
#' @param lo A non-negative integer. Minimum number of repeats, when grouped.
#' @param hi positive integer. Maximum number of repeats, when grouped.
#' @param lazy A logical value.  Should repetition be matched lazily or
#' greedily?
#' @param char_class A logical value. Should \code{x} be wrapped in a character
#' class?  If \code{NA}, the function guesses whether that's a good idea.
#' @return A character vector representing part or all of a regular expression.
#' @references \url{http://www.regular-expressions.info/repeat.html} and
#' \url{http://www.rexegg.com/regex-quantifiers.html}
#' @examples
#' # Can match constants or class values
#' repeated(GRAPH, 2, 5)
#' repeated(graph(), 2, 5)   # same
#'
#' # Short cuts for special cases
#' optional(blank())         # same as repeated(blank(), 0, 1)
#' zero_or_more(hex_digit()) # same as repeated(hex_digit(), 0, Inf)
#' one_or_more(printable())  # same as repeated(printable(), 1, Inf)
#'
#' # 'Lazy' matching (match smallest no. of chars)
#' repeated(cntrl(), 2, 5, lazy = TRUE)
#' lazy(one_or_more(cntrl()))
#'
#' # Overriding character class wrapping
#' repeated(ANY_CHAR, 2, 5, char_class = FALSE)
#' @include internal.R
#' @export
repeated <- function(x, lo, hi, lazy = FALSE, char_class = NA)
{
  x <- wrap_in_char_class(x, char_class)

  if(missing(lo))
  {
    lo <- 1L
  }
  lo <- as.integer(lo)
  if(missing(hi))
  {
    hi <- lo
  }
  args <- recycle(x, lo, hi, lazy)
  with(
    args,
    {
      hi <- ifelse(is.na(hi), lo, floor(hi))
      if(!all(is.finite(lo)))
      {
        stop("lo has missing or infinite values.")
      }
      if(any(lo < 0))
      {
        stop("lo has negative values.")
      }
      if(any(hi < lo))
      {
        stop("hi has values that are less than the corresponding values in lo.")
      }
      if(any(is.na(lazy)))
      {
        stop("lazy has missing values.")
      }
      rx <- regex(x, "{", lo, ",", hi, "}")

      # Special cases
      i <- lo == 0 & hi == 1
      rx[i] <- optional(x[i], char_class = FALSE)
      i <- lo == 1 & hi == 1
      rx[i] <- x[i]
      i <- lo == 0 & is.infinite(hi)
      rx[i] <- zero_or_more(x[i], char_class = FALSE)
      i <- lo == 1 & is.infinite(hi)
      rx[i] <- one_or_more(x[i], char_class = FALSE)
      i <- lo > 1 & is.infinite(hi)
      rx[i] <- regex(x[i], "{", lo[i], ",", "}")
      i <- lo != 1 & lo == hi
      rx[i] <- regex(x[i], "{", lo[i], "}")

      rx[lazy] <- lazy(rx[lazy])
      rx
    }
  )
}

#' @rdname repeated
#' @export
optional <- function(x, char_class = NA)
{
  x <- wrap_in_char_class(x, char_class)
  regex(x, "?")
}

#' @rdname repeated
#' @export
lazy <- function(x)
{
  regex(x, "?")
}

#' @rdname repeated
#' @export
zero_or_more <- function(x, char_class = NA)
{
  x <- wrap_in_char_class(x, char_class)
  regex(x, "*")
}

#' @rdname repeated
#' @export
one_or_more <- function(x, char_class = NA)
{
  x <- wrap_in_char_class(x, char_class)
  regex(x, "+")
}

wrap_in_char_class <- function(x, char_class = NA)
{
  if(is.na(char_class))
  {
    # Starts with '[' (not '[:'), and ends with ']' (not ':]')
    char_class <- !grepl("^\\[(?!:).*(?<!:)\\]$", x, perl = TRUE)
  }
  if(char_class)
  {
    x <- char_class(x)
  }
  x
}

