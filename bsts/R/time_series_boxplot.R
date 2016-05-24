# Copyright 2011 Google Inc. All Rights Reserved.
# Author: stevescott@google.com (Steve Scott)

TimeSeriesBoxplot <- function(x, time, ylim = NULL, add = FALSE, ...) {
  ## Plots side by side boxplots against horizontal times series axis.
  stopifnot(inherits(time, "Date") || inherits(time, "POSIXt"))
  if (is.null(ylim)) {
    ylim <- range(x)
  }
  if (!add) {
    plot(time, 1:length(time), type = "n", ylim = ylim)
  }
  dt <- mean(diff(time))
  boxplot(x, add = TRUE, at = time, show.names = FALSE,
          boxwex = dt/2, ylim = ylim, ...)
}
