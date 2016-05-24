# #########################################################################
# This file is a part of FCNN4R.
#
# Copyright (c) Grzegorz Klima 2015-2016
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
# #########################################################################


#' Plotting multilayer perceptron network
#'
#' This function plots a multilayer perceptron network's structure. Optionally,
#' weights' values are displayed on graph.
#'
#' @param net an object of \code{mlp_net} class
#' @param show_weights logical, should weights' values be displayed?
#'        (FALSE by default)
#' @param show_neuron_idx logical, should neurons' indices be displayed?
#'        (TRUE by default)
#'
#' @return This function does not return value.
#'
#' @export mlp_plot
#'
mlp_plot <- function(net, show_weights = FALSE, show_neuron_idx = TRUE)
{
    if (!is.mlp_net(net)) {
        stop("expected net argument to be of mlp_net class")
    }

    lays <- mlp_get_layers(net)
    L <- length(lays)
    dx <- max(max(lays) / L, 1)
    xl <- L * dx
    yl <- max(lays)
    layspts <- list()
    for (l in 1:L) {
        ly <- (1:lays[l])
        ly <- 0.5 * yl + ly - 0.5 * (1 + lays[l])
        layspts[[l]] <- rev(ly)
    }

    asp <- min(1, (yl / xl)^-.5)
    colio <- "black"
    colh <- "grey40"
    colw <- "blue4"

    plot.new()
    plot.window(xlim = c(0, xl), ylim = c(0, yl), asp = asp)

    ly <- layspts[[1]]
    for (n in 1:lays[1]) {
        x0 <- 0.5 * dx
        segments(x0 = x0 - .65, y0 = ly[n], x1 = x0 - .5, col = colio, lwd = 2)
        segments(x0 = x0 - .2, y0 = ly[n], x1 = x0, col = colio, lwd = 2)
        segments(x0 = x0 - .5, y0 = ly[n] - .05, x1 = x0 - .2, y1 = ly[n], col = colio, lwd = 2)
        segments(x0 = x0 - .5, y0 = ly[n] + .05, x1 = x0 - .2, y1 = ly[n], col = colio, lwd = 2)
        segments(x0 = x0 - .5, y0 = ly[n] + .05, y1 = ly[n] - .05, col = colio, lwd = 2)
    }
    ly <- layspts[[L]]
    for (n in 1:lays[L]) {
        x0 <- (l - .5) * dx
        segments(x0 = x0, y0 = ly[n], x1 = x0 + .3, col = colio, lwd = 2)
        segments(x0 = x0 + .3, y0 = ly[n] - .05, x1 = x0 + .6, y1 = ly[n], col = colio, lwd = 2)
        segments(x0 = x0 + .3, y0 = ly[n] + .05, x1 = x0 + .6, y1 = ly[n], col = colio, lwd = 2)
        segments(x0 = x0 + .3, y0 = ly[n] + .05, y1 = ly[n] - .05, col = colio, lwd = 2)
    }
    for (l in 2:L) {
        ly <- layspts[[l]]
        lyp <- layspts[[l - 1]]
        for (n in 1:lays[l]) {
            for (npl in 1:lays[l - 1]) {
                if (mlp_get_w_st(net, layer = l, nidx = n, nplidx = npl)) {
                    x0 <- (l - 1.5) * dx
                    y0 <- lyp[npl]
                    x1 <- (l - .5) * dx
                    y1 <- ly[n]
                    segments(x0 = x0, y0 = y0, x1 = x1, y1 = y1, col = colw)
                    if (show_weights) {
                        wv <- mlp_get_w(net, layer = l, nidx = n, nplidx = npl)
                        t <- .3
                        text(x = t * x0 + (1 - t) * x1, y = t * y0 + (1 - t) * y1,
                             labels = signif(wv, digits = 5), pos = 3, cex = .66,
                             srt = atan((y1 - y0) / dx * asp) / pi * 180, col = colw)
                    }
                }
            }
            if (mlp_get_w_st(net, layer = l, nidx = n, nplidx = 0)) {
                x0 <- (l - .5) * dx
                segments(x0 = x0, y0 = ly[n], y1 = ly[n] - .5, col = colw)
                    if (show_weights) {
                        wv <- mlp_get_w(net, layer = l, nidx = n, nplidx = 0)
                        text(x = x0, y = ly[n] - .1,
                             labels = signif(wv, digits = 5),
                             pos = 2, cex = .66, srt = 90, col = colw)
                }
            }
        }
    }
    for (l in 1:L) {
        ly <- layspts[[l]]
        if (l %in% c(1, L)) {
            cl <- colio
        } else {
            cl <- colh
        }
        points(x = rep((l - 0.5) * dx, length(ly)), y = ly, col = cl, pch = 16, cex = 5 / max(lays) ^ .3)
        if (show_neuron_idx) {
            text(x = rep((l - 0.5) * dx, length(ly)), y = ly, labels = 1:lays[l],
                cex = 1.3 / max(lays) ^ .3, col = "white")
        }
    }
}






