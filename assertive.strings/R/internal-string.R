# TODO: rework this using rebus

#' Create a regex from components.
#'
#' Creates a regex from regex components.
#' @param ... Character vectors of regex components.
#' @param l A list of character vectors for alternate specification.
#' @param sep Regex for separating components of complete regex.
#' Defaults to "an optional space or hyphen".
#' @return A string containing a regex.
#' Each element in the vectors are pasted together, separated by the
#' \code{sep} value.  Those character vectors are then preceded by "^"
#' (regex for 'start of string'() and followed by "$" (regex for end
#' of string).  Finally, the regexes are collapsed with "|" (regex for
#' 'or').
#' @examples
#' \dontrun{
#' cas_number_components <- c(
#'   "[[:digit:]]{1,7}", "[[:digit:]]{2}", "[[:digit:]]"
#' )
#' cas_number_rx <- create_regex(rx_components, sep = "-")
#' }
#' @importFrom assertive.base merge_dots_with_list
#' @importFrom assertive.base parenthesise
create_regex <- function (..., l = list(), sep = "[- ]?")
{
  x <- merge_dots_with_list(..., l = l)
  rx <- vapply(
    x,
    function(x)
    {
      parenthesise(paste0(parenthesise(x), collapse = sep))
    },
    character(1)
  )
  paste0("^", rx, "$", collapse = "|")
}

#' Create regex for repeated digits
#' 
#' Creates a regex string for repeated digits.
#' 
#' @param lo Minimum number of digits to match.
#' @param hi Optional maximum number of digits to match.
#' @param optional If \code{TRUE}, the digits are optional.
#' @note If \code{hi} is omitted, the returned regex will only match the exact 
#' number of digits given by \code{lo}.
#' @return A character vector of regexes.
#' @examples
#' \dontrun{
#' d(1:5)
#' d(1:5, 6:8)
#' d(0:2, Inf)
#' }
d <- function(lo, hi = NA_integer_, optional = FALSE)
{
  lo <- as.integer(lo)
  if(any(lo < 0))
  {
    stop("lo contains negative values.")
  }
  l <- recycle(lo = lo, hi = hi)
  lo <- l$lo
  hi <- l$hi
  rx <- ifelse(
    is.na(hi),
    {    
      sub("{1}", "", paste0("[[:digit:]]{", lo, "}"), fixed = TRUE)
    },
    {
      ifelse(
        is.infinite(hi) & hi > 0,
        {
          ifelse(
            lo == 0,
            {
              "[[:digit:]]*"
            },
            ifelse(
              lo == 1,
              {
                "[[:digit:]]+"
              },
              {
                paste0("[[:digit:]]{", lo, ",}")
              }
            )
          )
        },
        {
          hi <- as.integer(hi)
          if(any(hi <= lo))
          {
            stop("hi must be strictly greater than lo.")
          }
          rx <- paste0("[[:digit:]]{", lo, ",", hi, "}")
        }
      )
    } 
  )
  if(optional)
  {
    rx <- paste0("(", rx, ")?")
  }
  rx
}

#' Does the input match the regular expression?
#' 
#' Checks that the input matches the regular expression.
#'
#' @param x Input to check.
#' @param rx A regular expression.
#' @param ignore.case Should the case of alphabetic characters be ignored?
#' @param ... Passed to \code{\link{grepl}}.
#' @note The default for \code{ignore.case} is different to the default in 
#' \code{grepl}.
#' @return A logical vector that is \code{TRUE} when the input matches the 
#' regular expression.
#' @seealso \code{\link{regex}} and \code{\link{regexpr}}.
matches_regex <- function(x, rx, ignore.case = TRUE, ...)
{
  call_and_name(
    function(x) 
    {
      if(!nzchar(rx[1]))
      {
        warning(
          "Regular expression is the empty string, and matches everything."
        )
        return(rep.int(TRUE, length(x)))
      }
      # call to ifelse needed because grepl always returns TRUE or FALSE
      # need to unname, because ifelse preserves x's names, when we want to
      # name result with values of x, and merge.list throws a warning about
      # duplicate names attr.
      ifelse(   
        is.na(unname(x)),
        NA,
        grepl(rx, x, ignore.case = ignore.case, ...)
      )
    }, 
    x
  )
}

#' Recycle arguments
#' 
#' Explicit recycling of arguments to make them all have the same length.
#' @param ... Arguments, usually vectors.
#' @return A \code{list} of vectors, all with the same length.
#' @note The function is based on \code{rep_len}, which drops attributes (hence
#' this being most appropriate for vector inputs).
#' @seealso \code{\link[base]{rep_len}}.
#' @examples
#' \dontrun{
#' # z is the longest argument, with 6 elements
#' recycle(x = 1:4, y = list(a = month.abb, b = pi), z = matrix(1:6, nrow = 3))
#' }
recycle <- function(...)
{
  dots <- list(...)
  n <- max(vapply(dots, length, integer(1)))
  lapply(dots, rep_len, length.out = n)
}

#' Removes invalid characters from a string.
#'
#' Removes invalid characters from a string, leaving only digits.
#' @param x Input to strip.
#' @param invalid_chars A regular expression detailing characters to remove.
#' @param char_desc A string describing the characters to remove.
#' @param allow_x If \code{TRUE}, the letter "X" is allowed - useful for check 
#' digits.
#' @param allow_plus If \code{TRUE}, the symbol "+" is allowed - useful for 
#' phone numbers.
#' @return A character vector of the same length as \code{x}, consisting of 
#' strings without the characters detailed in the \code{invalid_chars}.
#' @examples
#' \dontrun{
#' strip_invalid_chars(
#'   "  We're floating\tin    space\n\n\n", "[[:space:]]", "whitespace"
#' )
#' strip_non_numeric(" +44 800-123-456 ", allow_plus = TRUE)
#' #Inputs such as factors as coerced to character.
#' strip_non_alphanumeric(factor(c(" A1\t1AA.", "*(B2^2BB)%")))
#' }
strip_invalid_chars <- function(x, invalid_chars, char_desc = gettext("invalid"))
{
  x <- coerce_to(x, "character")
  if(any(grepl(invalid_chars, x)))
  {
    warning(gettextf("Removing %s characters from input.", char_desc))
    x <- gsub(invalid_chars, "", x)
  }
  x
}

#' @rdname strip_invalid_chars
strip_non_alphanumeric <- function(x)
{
  strip_invalid_chars(x, "[^[:alnum:]]+", "non-alphanumeric")
}

#' @rdname strip_invalid_chars
strip_non_numeric <- function(x, allow_x = FALSE, allow_plus = FALSE)
{
  invalid_chars <- paste0(
    "[^[:digit:]", 
    if(allow_x) "X", 
    if(allow_plus) "\\+", 
    "]+", 
    collapse = ""
  )
  strip_invalid_chars(x, invalid_chars, "non-numeric")
}
