#' Generate random fake latitude values.
#'
#' @param n number of observations.
#' @return A character vector of \code{n} randomly generated latitude values.
#' @examples
#' r_latitudes(10)
#' @export
r_latitudes <- function(n) {
  stats::runif(n, -90, 90)
}
