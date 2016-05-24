#' Escape special characters
#' 
#' Prefix the special characters with a blackslash to make them literal 
#' characters.
#' @param x A character vector.
#' @param escape_brace A logical value indicating if opening braces should be 
#' escaped.  If using R's internal PRCE engine or \code{stringi}'s ICU engine, 
#' you want this.  If using the perl engine, you don't.
#' @return A character vector, with regex meta-characters escaped.
#' @note Special characters inside character classes (except caret, hypen and 
#' closing bracket in certain positions) do not need to be escaped. This 
#' function makes no attempt to parse your regular expression and decide whether 
#' or not the special character is inside a character class or not.  It simply
#' escapes every value.
#' @examples
#' escape_special("\\ ^ $ . | ? * + ( ) { } [ ]")
#' @export
escape_special <- function(x, escape_brace = TRUE)
{
  specials <- c("\\", "^", "$", ".", "|", "?", "*", "+", "(", ")", "[", if(escape_brace) "{")
  for(char in specials)
  {
    x <- gsub(char, paste0("\\", char), x, fixed = TRUE)
  }
  regex(x)
}
