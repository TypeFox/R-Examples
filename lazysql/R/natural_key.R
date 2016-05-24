#' Create SQL string for joining on matching natural keys
#'
#' @description
#'  Create string with SQL expressions for \code{WHERE} clause
#'  to join two tables on the given columns.
#'
#' @param table_names [character(2)]\cr
#'  Name of data base tables to be joined.
#' @param key_columns [character(1:Inf)]\cr
#'  Names of key columns in both tables.
#'
#' @details
#'  The names of tables and key columns must be valid SQL identifiers.
#'  They are validated to conform to
#'  the regular expression returned by \code{\link{valid_identifier_regex}}.
#'
#'  The SQL string is created in 3 steps:
#'  \enumerate{
#'    \item Combine table names with key names, eg, "\code{PRL.FLIGHT_NR}".
#'    \item Create logical expressions, eg, "\code{PRL.FLIGHT_NR = PRL_SSR.FLIGHT_NR}"
#'    \item Concatenate logical expressions by \code{"and"} to form final SQL esxpression.
#'  }
#'
#' @note
#'  The current implementation assumes that key columns have the same names in both tables.
#'
#' @return
#'  Character string to be used in SQL statement.
#' @author Uwe Block
#'
#' @examples
#' # SQL expression
#' (sql_expr <- lazysql::natural_key(c("TAB1", "tab_2"),c("COL1", "col_2")))
#'
#' # sample SQL JOIN statement
#' paste("select * from TAB1, TAB2 where", sql_expr)
#'
#' @seealso \code{\link{valid_identifier_regex}}.
#' @import magrittr
#' @export
natural_key <- function(
  table_names,
  key_columns
) {
  checkmate::assert_character(table_names, min.chars = 1L, len = 2L, any.missing = FALSE,
                              unique = TRUE, pattern = valid_identifier_regex())
  checkmate::assert_character(key_columns, min.chars = 1L, min.len = 1L, any.missing = FALSE,
                              unique = TRUE, pattern = valid_identifier_regex())
  sql <-
    vapply(table_names, FUN = paste, FUN.VALUE = character(length(key_columns)),
           key_columns, sep = ".") %>%
    plyr::aaply(.margins = 1, .fun = paste, collapse = " = ") %>%
    paste(collapse = " and ")
  return(sql)
}
