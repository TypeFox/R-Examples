#' Create SQL string to select date between two given dates
#'
#' @description
#'  Create string with SQL \code{BETWEEN} expression for \code{WHERE} clause to select dates
#'  within the given range.
#'
#' @param column_name [character(1)]\cr
#'  Name of data base column to select dates from.
#' @param date_range [Date(1:2)]\cr
#'  One or two dates giving the date range in which the dates should be enclosed (closed interval).
#'  If only one date is given, it is taken for both upper and lower limits.

#' @details
#'  \code{column_name} must be a valid SQL identifier. It is validated to conform to
#'  the regular expression returned by \code{\link{valid_identifier_regex}}.
#'
#' @return
#'  Character string to be used in SQL statement.

#' @author Uwe Block

#' @examples
#' date1 <- as.Date("2016-02-22")
#' date2 <- as.Date("2016-02-11")
#'
#' # SQL expression for a date range
#' (sql_expr1 <- lazysql::date_between("STD_1", c(date1, date2)))
#'
#' # SQL expression for a single date
#' (sql_expr2 <- lazysql::date_between("STD_1", date1))
#'
#' # sample SQL statements
#' paste("select * from TEST_TABLE where", sql_expr1)
#'
#' paste("select * from TEST_TABLE where", sql_expr2)
#'
#' @seealso \code{\link{valid_identifier_regex}}.
#' @export
date_between <- function(
  column_name,
  date_range
) {
  checkmate::assert_string(column_name, na.ok = FALSE, min.chars = 1L,
                           pattern = valid_identifier_regex())
  checkmate::assert_numeric(date_range, any.missing = FALSE,
                            min.len = 1L, max.len = 2L)
  checkmate::assert_class(date_range, "Date")
  fmt <- "to_date('%s', 'yyyy-mm-dd')"
  date_from <- sprintf(fmt, min(date_range))
  date_to <- sprintf(fmt, max(date_range))
  sql <- paste(column_name, "between", date_from, "and", date_to)
  return(sql)
}
