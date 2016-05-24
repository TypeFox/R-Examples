#' Creates the limits
#' \code{limits} creates the limits
#' @keywords internal
#' @export
limits <- function(d, x, xmin, xmax, log) {
  d %>% do(one_limit(., x, xmin, xmax, log))
}
