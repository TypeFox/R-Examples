#' Creates the limits for one condition
#' \code{one_limit} creates the limits for one condition
#' @keywords internal
#' @export
one_limit <- function(d, x, xmin, xmax, log) {
  if (is.null(xmin)) xmin <- min(d[[x]])
  if (is.null(xmax)) xmax <- max(d[[x]])
  data.frame(xmin, xmax)
}
