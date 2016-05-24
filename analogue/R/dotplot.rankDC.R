dotplot.rankDC <- function(x, data = NULL, sort = TRUE, decreasing = FALSE,
                           xlab = "Rank correlation", ...) {
  if(sort)
    x <- sort(x, decreasing = decreasing)
  dotplot(~ x, xlab = xlab, ...)
}
