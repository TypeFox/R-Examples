DiffMeans <- function (y, tr, tl, ...) {
  ## Difference of means
  mean.default(y[tr == tl]) - mean.default(y[tr != tl])
}
