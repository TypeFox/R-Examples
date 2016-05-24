#
# file tree/R/plot.tree.sequence.R copyright (C) 1994-2012 B. D. Ripley
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#
plot.tree.sequence <- function(x, ..., type = "l", ylim = range(x$dev),
    order = c("increasing", "decreasing"))
{
    if(missing(type) && inherits(x, "prune")) type <- "S"
    if(is.null(x$method)) x$method <- "deviance"
    order <-  match.arg(order)
    if(order == "increasing")
        sign <- +1 else sign <- -1
    plot(sign*x$size, x$dev, axes = FALSE,
         xlab = "size", ylab = x$method,
         type = type, ylim = ylim, ...)
    box()
    axis(2L, ...)
    xaxp <- par("xaxp")
    pos <- sign*seq(xaxp[1L], xaxp[2L], diff(xaxp[-3])/xaxp[3L])
    if(pos[1L] == 0) pos[1L] <- 1
    n <- length(pos)
    maxsize <- max(x$size)
    if(pos[n] > maxsize) pos[n] <- maxsize
    axis(1L, at = sign*pos, labels = pos, ...)
    axis(3L, at = sign * x$size, labels = format(signif(x$k, 2L)), ...)
    invisible()
}
