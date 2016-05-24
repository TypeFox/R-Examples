#' Generate random fake date of birth values.
#'
#' @param n number of observations.
#' @param start starting date.
#' @param end ending date.
#' @return A character vector of \code{n} randomly generated date of birth values.
#' @examples
#' r_date_of_births(10)
#' r_date_of_births(10, start = as.Date("2000-01-01"))
#' r_date_of_births(10,
#'                  start = as.Date("2000-01-01"),
#'                  end = as.Date("2100-01-01"))
#' @export
r_date_of_births <- function(n, start = as.Date("1900-01-01"), end = Sys.Date()) {
  as.Date(sample.int(end - start, size = n), origin = start)
}
