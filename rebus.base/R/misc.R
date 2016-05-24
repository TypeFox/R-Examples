#' Make a regex exact
#'
#' Makes a regex exact: that is, it must contain the whole string, not just part
#' of it.
#' @param x A character vector.
#' @return A character vector representing part or all of a regular expression.
#' @examples
#' # A hex color
#' (rx <- "#" %R% hex_digit(6))
#' (rx_exact <- exactly(rx))
#'
#' # Usage
#' stringi::stri_detect_regex("ginger is #B06500", rx)
#' stringi::stri_detect_regex("ginger is #B06500", rx_exact)
#' stringi::stri_detect_regex("#B06500", rx_exact)
#' @export
exactly <- function(x)
{
  regex(START, x, END)
}

#' Treat part of a regular expression literally.
#'
#' Treats its contents as literal characters.
#' @param x A character vector.
#' @return A character vector representing part or all of a regular expression.
#' @examples
#' (rx <- digit(1, 3))
#' (rx_literal <- literal(rx))
#'
#' # Usage
#' stringi::stri_detect_regex("123", rx)
#' stringi::stri_detect_regex("123", rx_literal)
#' stringi::stri_detect_regex("[[:digit:]]{1,3}", rx_literal)
#' @export
literal <- function(x)
{
  regex("\\Q", x, "\\E")
}

