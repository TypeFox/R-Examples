#' Convert All Factor Columns to Character Columns
#'
#' By default, R converts character columns to factors.
#' Instead of re-reading the data using \code{stringsAsFactors}, the
#' \code{\link{safe.chars}} function will identify which columns are currently factors, and convert them all to characters.
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}.
#'
#' @param .data The name of the \code{data.frame}
#' @seealso \code{\link{read.table}}, \code{\link{destring}}.
#' @examples
#'  str(iris)
#' iris_2 = safe.chars(iris)
#' str(iris_2)
#'
#' @export
safe.chars <- function(.data) {
  .data[sapply(.data, is.factor)] <-
    lapply(.data[sapply(.data, is.factor)], as.character)
  .data
}
