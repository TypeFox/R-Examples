#' Convert Proportions to the Probit Scale
#'
#' Convert proportions to the probit scale.
#' @param prob
#'   A numeric vector of proportions.
#' @return
#'   A numeric vector the same length as \code{prob} with quantiles on
#'     the probit scale.
#' @export
#' @details
#'   Simply calls \code{\link{qnorm}(prob)}.
#' @examples
#' probit(c(0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999))

probit <- function(prob) {
  qnorm(prob)
}
