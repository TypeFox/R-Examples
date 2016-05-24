#' Generate random fake phone numbers.
#'
#' @param n number of observations.
#' @param use_hyphens should hyphens be included.
#' @param use_parentheses should parantheses be included.
#' @param use_spaces should spaces be included.
#' @return A character vector of \code{n} randomly generated phone numbers.
#' @examples
#' r_phone_numbers(10)
#' r_phone_numbers(10, use_hyphens = TRUE)
#' r_phone_numbers(10, use_parentheses = TRUE)
#' r_phone_numbers(10, use_spaces = TRUE)
#' r_phone_numbers(10, use_parentheses = TRUE, use_hyphens = TRUE)
#' r_phone_numbers(10, use_parentheses = TRUE, use_spaces = TRUE)
#' @export
r_phone_numbers <- function(n, use_hyphens = FALSE, use_parentheses = FALSE, use_spaces = FALSE) {
  left_paren <- ""
  right_paren <- ""
  if(use_parentheses) {
    left_paren <- "("
    right_paren <- ")"
  }
  hyphen <- ""
  if(use_hyphens) hyphen <- "-"
  space <- ""
  if(use_spaces) space <- " "
  build_phone_number <- function(l, r, h, s) {
    paste0(l, paste0(sample(1:9, size = 3), collapse = ""), r, h, s,
           paste0(sample(1:9, size = 3), collapse = ""), h, s,
           paste0(sample(1:9, size = 4), collapse = ""),
           collapse = "")
  }


  return(replicate(n, build_phone_number(l = left_paren,
                                         r = right_paren,
                                         h = hyphen,
                                         s = space)))
}
