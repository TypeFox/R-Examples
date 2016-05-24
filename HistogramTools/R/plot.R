# Copyright 2013 Google Inc. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Author: mstokely@google.com (Murray Stokely)

PlotRelativeFrequency <- function(x, ylab="Relative Frequency", ...) {
# Plots relative frequency histogram.
#
# By default R plots only absolute frequency or density.
#
# Args:
#   x:    A histogram
#   ...:  Additional arguments to pass to plot()
  x <- ScaleHistogram(x, 1 / Count(x))
  plot(x, freq=TRUE, ylab=ylab, ...)
}

PlotLog2ByteEcdf <- function(x,
                             xlab="Bytes (log)",
                             ylab="Cumulative Fraction",
                             with.grid=TRUE,
                             ...) {
# Plots ECDF of a histogram with power of 2 byte boundaries.
#
# Args:
#   x: A histogram or ecdf object.
#   with.grid: If TRUE, draw a faint grid on the plot.
#   cex.axis: cex parameter for the axes.
#   ...: Additional arguments to pass to plot()
  if (inherits(x, "histogram")) {
    x <- HistToEcdf(x)
  }
  stopifnot(inherits(x, "ecdf"))
  power.of.two.breaks <- intersect(knots(x), 2^(0:53))
  if (length(power.of.two.breaks) < 3) {
    stop("Insufficient powers of 2 in knots() of ecdf")
  }
  plot(knots(x), x(knots(x)),
       type="l",
       lwd=2,
       log="x",
       xaxt="n",
       las=1,
       xlab=xlab,
       ylab=ylab,
       ...
       )
  if (requireNamespace("gdata", quietly=TRUE)) {
    labs <- gdata::humanReadable(power.of.two.breaks)
  } else {
    labs <- power.of.two.breaks
  }
  axis(1, at=power.of.two.breaks, labels=labs)
  if (with.grid) {
    abline(h=seq(.2,.8, by=.2), lty="dotted", col="lightgray")
    abline(v=power.of.two.breaks, lty="dotted", col="lightgray")
  }
}

PlotLogTimeDurationEcdf <- function(x, with.grid=TRUE,
                                    xlab="Age (log)",
                                    ylab="Cumulative Fraction",
                                    cex.lab=1.6,
                                    cex.axis=1.6, ...) {
# Plots weighted ECDF of a histogram with log time bucket boundaries
# (seconds since the epoch) with log-scaled x-axis labels from 1 minute
# to 1 year.
#
# Args:
#   x: A histogram or ecdf object.
#   with.grid: If TRUE, draw a faint grid on the plot.
#   cex.axis: cex parameter for the axes.
#   ...: Additional arguments to pass to plot()

  if (inherits(x, "histogram")) {
    x <- HistToEcdf(x)
  }
  stopifnot(inherits(x, "ecdf"))

  plot(knots(x), x(knots(x)),
       type="l",
       log="x",
       xaxt="n",
       xlab=xlab,
       ylab=ylab,
       yaxt="n",
       cex.lab=cex.lab,
       ...)

  xticks <- list("1s" = 1,
                 "10s" = 10,
                 "1m" = 60,
                 "10m" = 60 * 10,
                 "1h" = 3600,
                 "6h" = 3600 * 6,
                 "1d" = 86400,
                 "3d" = 86400 * 3,
                 "7d" =  86400 * 7,
                 "30d" = 86400 * 30,
                 "90d" = 86400 * 30 * 3,
                 "1y" = 86400 * 365,
                 "10y" = 86400 * 365 * 10)

  axis(1, at=unlist(xticks), labels=names(xticks), cex.axis=cex.axis)
  axis(2, at=seq(0, 1, by=.2), las=1, cex.axis=cex.axis)

  if (with.grid) {
    abline(h=seq(.2,.8, by=.2), lty="dotted", col="lightgray")
    abline(v=xticks, lty="dotted", col="lightgray")
  }
}
