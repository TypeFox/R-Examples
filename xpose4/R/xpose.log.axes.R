# Xpose 4
# An R-based population pharmacokinetic/
# pharmacodynamic model building aid for NONMEM.
# Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
# Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
# Mats O. Karlsson and E. Niclas Jonsson.
# Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
# E. Niclas Jonsson.

# This file is a part of Xpose 4.
# Xpose 4 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  A copy can be cound in the R installation
# directory under \share\licenses. If not, see http://www.gnu.org/licenses/.

xpose.logTicks <- function (lim, loc = c(1, 5)) {
  ii <- floor(log10(range(lim))) + c(-1, 2)
  main <- 10^(ii[1]:ii[2])
  r <- as.numeric(outer(loc, main, "*"))
  r[lim[1] <= r & r <= lim[2]]
}  

xpose.yscale.components.log10 <- function(lim, ...) {
  ans <- yscale.components.default(lim = lim, ...)
  tick.at <- xpose.logTicks(10^lim, loc = 1:9)
  tick.at.major <- xpose.logTicks(10^lim, loc = 1)
  major <- tick.at %in% tick.at.major
  ans$left$ticks$at <- log(tick.at, 10)
  ans$left$ticks$tck <- ifelse(major, 1.5, 0.75)
  ans$left$labels$at <- log(tick.at, 10)
  ans$left$labels$labels <- as.character(tick.at)
  ans$left$labels$labels[!major] <- ""
  ans$left$labels$check.overlap <- FALSE
  ans
}

xpose.xscale.components.log10 <- function(lim, ...) {
  ans <- xscale.components.default(lim = lim, ...)
  tick.at <- xpose.logTicks(10^lim, loc = 1:9)
  tick.at.major <- xpose.logTicks(10^lim, loc = 1)
  major <- tick.at %in% tick.at.major
  ans$bottom$ticks$at <- log(tick.at, 10)
  ans$bottom$ticks$tck <- ifelse(major, 1.5, 0.75)
  ans$bottom$labels$at <- log(tick.at, 10)
  ans$bottom$labels$labels <- as.character(tick.at)
  ans$bottom$labels$labels[!major] <- ""
  ans$bottom$labels$check.overlap <- FALSE
  ans
}
