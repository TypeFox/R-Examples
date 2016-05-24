#' Generate random fake longitude values.
#'
#' @param n number of observations.
#' @return A character vector of \code{n} randomly generated longitude values.
#' @examples
#' r_longitudes(10)
#' @export
r_longitudes <- function(n) {
  stats::runif(n, -180, 180)
}
