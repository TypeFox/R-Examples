# avoid CRAN notes when using magrittr
utils::globalVariables(".")

#' Regex pattern to validate SQL identifier names
#'
#' @description
#' Returns a regular expression to validate unquoted SQL identifiers.
#' @details
#' Valid SQL identifiers must begin with an alphabetic character followed by
#' alphanumeric characters or underscores "\code{_}".
#' @note
#' The current implementation doesn't allow any other special characters in
#' SQL identfiers or quoted SQL identifiers for safety reasons.
#' In future releases, valid SQL identifiers might be defined depending
#' on the target database system.
#' @return Character string with regular expression.
#' @references
#' ORACLE Database SQL Language Reference.
#' @author Uwe Block
#' @examples
#' lazysql::valid_identifier_regex()
#' @export
#'
valid_identifier_regex <- function() {
  return("^[[:alpha:]][_[:alnum:]]*$")
}
