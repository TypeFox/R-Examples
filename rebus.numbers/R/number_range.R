#' @include imports.R
NULL

#' Generate a regular expression for a number range.
#'
#' Generates a regular expression that matches a sequence of numbers.
#' @param lo An integer.
#' @param hi An integer greater than or equal to \code{lo}.
#' @param allow_leading_zeroes A logical value. Are leading zeroes allowed to
#' bulk the match up to the length of the number with the most digits?
#' @param capture A logical value. See \code{\link{or}} for details.
#' @return A character vector representing part or all of a regular expression.
#' @examples
#' number_range(0, 255)
#' number_range(0, 255, allow_leading_zeroes = TRUE)
#' number_range(10000, 19999)
#' number_range(6, 54321)
#' number_range(-77, 77)
#' number_range(-77, 77, capture = TRUE)
#' number_range(-77, 77, capture = NA)
#' @export
number_range <- function(lo, hi, allow_leading_zeroes = FALSE, capture = FALSE)
{
  lo <- as.integer(lo[1])
  hi <- as.integer(hi[1])
  if(hi < lo)
  {
    warning("hi is greater than lo.")
    return(UNMATCHABLE)
  }
  if(hi == lo)
  {
    return(lo)
  }
  if(lo < 0L && hi < 0L)
  {
    return("-" %R% number_range(-hi, -lo, allow_leading_zeroes, capture))
  }
  if(lo < 0L && hi > 0L)
  {
    return(
      engroup(
        number_range(lo, -1, allow_leading_zeroes) %|%
          number_range(0, hi, allow_leading_zeroes),
        capture = capture
      )
    )
  }
  x <- format(seq.int(lo, hi, 1L))
  d <- as.data.frame(
    do.call(rbind, strsplit(x, " ?")),
    stringsAsFactors = FALSE
  )
  alternates <- get_alternate_ranges(d, allow_leading_zeroes)
  rx <- simplify_repeated_digits(
    simplify_leading_zeroes(
      or1(alternates, capture = capture)
    )
  )
  as.regex(rx)
}

get_alternate_ranges <- function(d, allow_leading_zeroes)
{
  if(is.null(ncol(d)) || ncol(d) == 0)
  {
    return(NULL)
  }
  if(ncol(d) == 1)
  {
    return(char_range(d[1, 1], d[nrow(d), 1]))
  }
  if(max(d[, 1]) == min(d[, 1]))
  {
    return(
      max(d[, 1]) %R%
        get_alternate_ranges(d[, -1, drop = FALSE], allow_leading_zeroes)
    )
  }
  grp <- factor(
    ifelse(
      d[, 1] == min(d[, 1]) & (!all(d[1, -1] %in% c("", "0"))),
      "min",
      ifelse(
        d[, 1] == max(d[, 1]) & any(d[nrow(d), -1] != "9"),
        "max",
        "middle"
      )
    ),
    levels = c("min", "middle", "max")
  )
  split_d <- split(d, grp)
  with(
    split_d,
    c(
      if(nrow(min) > 0)
      {
        prefix <- if(min[1, 1] == "" && allow_leading_zeroes)
        {
          optional(0)
        } else
        {
          min[1, 1]
        }
        prefix %R%
          get_alternate_ranges(min[, -1, drop = FALSE], allow_leading_zeroes)
      } else
      {
        NULL
      },
      if(nrow(middle) > 0)
      {
        m1 <- if(middle[1, 1] == "") "0" else middle[1, 1]
        m2 <- if(middle[nrow(middle), 1] == "") "0" else middle[nrow(middle), 1]
        suppressWarnings(char_range(m1, m2)) %R%
          ascii_digit(ncol(d) - 1, ncol(d) - 1)
      } else
      {
        NULL
      },
      if(nrow(max) > 0)
      {
        max[1, 1] %R%
          get_alternate_ranges(max[, -1, drop = FALSE], allow_leading_zeroes)
      } else
      {
        NULL
      }
    )
  )
}

# These next two functions are slow.  Maybe speed up with trickery using strsplit?
simplify_repeated_digits <- function(x)
{
  if(length(x) > 1)
  {
    warning("Only using the first element of x.")
    x <- x[1]
  }
  rx <- "(\\Q[0-9]\\E){2,}"
  repeat
  {
    m <- regexpr(rx, x)
    if(m == -1) break
    match_len <- attr(m, "match.length")
    n <- match_len / 5
    x <- paste0(substring(x, 1, m - 1), ascii_digit(n), substring(x, m + match_len))
  }
  rx <- "\\Q[0-9]\\E(\\{[0-9]\\}){2}"
  repeat
  {
    m <- regexpr(rx, x)
    if(m == -1) break
    match_len <- attr(m, "match.length")
    n <- (match_len - 5) / 3
    x <- paste0(
      substring(x, 1, m - 1),
      "[0-9]{",
      as.integer(substring(x, m + 6, m + 6)) + as.integer(substring(x, m + 9, m + 9)) - 1,
      "}",
      substring(x, m + match_len)
    )
  }
  x
}

# @examples
# x <- "x[0-9]x[0-9][0-9]x[0-9][0-9][0-9]x0?x0?0?x0?0?0?x"
# simplify_leading_zeroes(x)
simplify_leading_zeroes <- function(x)
{
  if(length(x) > 1)
  {
    warning("Only using the first element of x.")
    x <- x[1]
  }
  rx <- "(0\\?){2,}"
  repeat
  {
    m <- regexpr(rx, x)
    if(m == -1) break
    match_len <- attr(m, "match.length")
    n <- match_len / 2
    x <- paste0(substring(x, 1, m - 1), repeated("0", 0, n), substring(x, m + match_len))
  }
  x
}
