#' @title Weighted Variance
#'
#' @description Weighted Variance Formula
#'
#' @param x the variable.
#' @param w the variance.
#' @param na.rm A logical if NA should be disregarded.
#' @keywords Stats
#' @export
#' @examples
#' wt=c(1.23, 2.12, 1.23, 0.32, 1.53, 0.59, 0.94, 0.94, 0.84, 0.73)
#' x = c(5, 5, 4, 4, 3, 4, 3, 2, 2, 1)
#' wtd.var(x, wt)
`wtd.var` <- function(x, w, na.rm = FALSE) {
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  sum.w <- sum(w)
  sum.w2 <- sum(w^2)
  mean.w <- sum(x * w) / sum(w)
  (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm = na.rm)
}
NULL
