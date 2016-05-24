#' Calculates the sum of squared errors of prediction
#' \code{one_sse} calculates  the sum of squared errors of prediction for one
#' condition
#' @keywords internal
#' @export
one_sse <- function(d, groups, averages) {
  if (length(groups) != 0) averages <- semi_join(averages, d, by = groups)
  data.frame(sse = sum((averages$y-d$ypred)^2))
}



