#' Create SQL string to select values included in a set of given values
#'
#' @description
#'  Create string with SQL \code{IN} expression for \code{WHERE} clause to select values
#'  included in a set of given values.
#'
#' @param column_name [character(1)]\cr
#'  Name of data base column to select values from.
#' @param choices [character(1:Inf)] or [integer(1:Inf)]\cr
#'  The values which must be matched. Character values must not contain any
#'  single or double quotes to avoid problems with SQL syntax and for safety reasons.
#' @param negation [character(1)]\cr
#'  If \code{"not"} the selection is inverted to a \code{NOT IN} expression.
#'
#' @details
#'  \code{column_name} must be a valid SQL identifier. It is validated to conform to
#'  the regular expression returned by \code{\link{valid_identifier_regex}}.
#'
#' @return
#'  Character string to be used in SQL statement.
#'
#' @author Uwe Block
#'
#' @examples
#' # SQL expressions
#' lazysql::in_condition("COL_1", 1:3)
#'
#' lazysql::in_condition("COL_1", 1:3, "not")
#'
#' lazysql::in_condition("COL_1", LETTERS[2:3])
#'
#' lazysql::in_condition("COL_1", LETTERS[2:3], "not")
#'
#' @seealso \code{\link{valid_identifier_regex}}.
#' @import magrittr
#' @export
in_condition <- function(
  column_name,
  choices,
  negation = c("", "not")
) {
  checkmate::assert_string(column_name, na.ok = FALSE, min.chars = 1L,
                           pattern = valid_identifier_regex())
  checkmate::assert(
    checkmate::checkCharacter(choices, any.missing = FALSE,
                              min.len = 1L, pattern = "^[^'\"]*$"),
    # see http://stackoverflow.com/a/198810 for the regex pattern
    checkmate::checkInteger(choices, any.missing = FALSE,
                            min.len = 1L)
  )
  checkmate::assert_subset(negation, eval(formals()$negation))
  negation <- match.arg(negation)
  checkmate::assert_string(negation, na.ok = FALSE, min.chars = 0L)
  # helper function
  prepare_values <- function(x) {
    if (is.character(x)) {
      return(paste0("'", x, "'"))
    } else {
      return(as.character(x))
    }
  }
  # build sql
  sql <-
    choices %>%
    prepare_values %>%
    paste(collapse = ", ") %>%
    paste0("(", ., ")") %>%
    paste(column_name, negation, "in", .)
  return(sql)
}
