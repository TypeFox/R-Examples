#' Apply mode modifiers
#' 
#' Applies one or more mode modifiers to the regular expression.
#' @param x A character vector.
#' @param modes A character vector of mode modifiers.
#' @return A character vector representing part or all of a regular expression.
#' @references \url{http://www.regular-expressions.info/modifiers.html} and
#' \url{http://www.rexegg.com/regex-modifiers.html}
#' @examples
#' x <- "foo"
#' case_insensitive(x)
#' free_spacing(x)
#' single_line(x)
#' multi_line(x)
#' duplicate_group_names(x)
#' no_backslash_escaping(x)
#' modify_mode(x, c("i", "J", "X"))
#' @include regex-methods.R
#' @export
modify_mode <- function(x, modes = c("i", "x", "s", "m", "J", "X"))
{
  modes <- paste0(match.arg(modes, several.ok = TRUE), collapse = "")
  regex("(?", modes, ")", x, "(?-", modes, ")") 
}

#' @rdname modify_mode
#' @export
case_insensitive <- function(x)
{
  modify_mode(x, "i")
}

#' @rdname modify_mode
#' @export
free_spacing <- function(x)
{
  modify_mode(x, "x")
}

#' @rdname modify_mode
#' @export
single_line <- function(x)
{
  modify_mode(x, "s")
}

#' @rdname modify_mode
#' @export
multi_line <- function(x)
{
  modify_mode(x, "m")
}

#' @rdname modify_mode
#' @export
duplicate_group_names <- function(x)
{
  modify_mode(x, "J")
}

#' @rdname modify_mode
#' @export
no_backslash_escaping <- function(x)
{
  modify_mode(x, "X")
}
