# This is a hidden function of the l2boost package.

# nice standard errors for plots
#
# @param  x Vector of error bar x value locations
# @param  upper Vector of upper error bar limits
# @param  lower Vector of limit error bar limits
# @param  width errorbar line width (default: 0.001)
# @param  max.M maximum number of bars to show in a plot (default: 100)
# @param  ... Additional arguments passed to segment function
#
# @seealso \code{\link{segments}}
#
error.bars <- function(x, upper, lower, width = 0.001, max.M = 100, ...) {
  M <- length(x)
  # thin-out x for presentable plots
  if (M > max.M) {
    pt <- unique(round(seq.int(1, max(x), length.out = max.M)))
  }
  else {
    pt <- 1:M
  }
  xlim <- range(x[pt])
  barw <- diff(xlim) * width
  segments(x[pt], upper[pt], x[pt], lower[pt], ...)
  segments(x[pt] - barw, upper[pt], x[pt] + barw, upper[pt], ...)
  segments(x[pt] - barw, lower[pt], x[pt] + barw, lower[pt], ...)
}
