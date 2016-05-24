#' Sum of squared errors of prediction
#'
#' \code{ypred} calculates the sum of squared errors of prediction
#' @param qp output from quickpsy
#' @export
sse <- function(qp) {
  qp$ypred %>% do(one_sse(., qp$groups, qp$averages))
}
