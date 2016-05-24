# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2    
"plot.surface" <- function(x, main = NULL, type = "C", 
    zlab = NULL, xlab = NULL, ylab = NULL, levels = NULL, zlim = NULL, 
    graphics.reset = NULL, labcex = 0.6, add.legend = TRUE, ...) {
    obj <- x
    old.par <- par(no.readonly = TRUE)
    if (is.na(match(type, c("b", "c", "C", "I", "p")))) {
        stop("plot type does not match b, C, I, or p.")
    }
    if (is.null(zlim)) {
        zlim = range(obj$z, na.rm = TRUE)
    }
    if (is.null(graphics.reset) & (type == "b")) {
        graphics.reset <- TRUE
    }
    else {
        graphics.reset <- FALSE
    }
    if (graphics.reset) {
        on.exit(par(old.par))
    }
    if (is.null(xlab)) {
        if (is.null(obj$xlab)) 
            xlab <- "X"
        else xlab <- obj$xlab
    }
    if (is.null(ylab)) {
        if (is.null(obj$ylab)) 
            ylab <- "Y"
        else ylab <- obj$ylab
    }
    if (is.null(zlab)) {
        if (is.null(obj$zlab)) 
            zlab <- "Z"
        else zlab <- obj$zlab
    }
    if (is.null(main)) 
        if (!is.null(obj$main)) 
            main <- obj$main
    if (type == "b") 
        set.panel(1, 2, TRUE)
    if (type == "p" | type == "b") {
        if (type == "b") {
            add.legend <- FALSE
            old.mar <- par()$mar
            par(mar = c(0, 5, 0, 0))
        }
        drape.plot(obj, xlab = xlab, ylab = ylab, zlab = zlab, 
            zlim = zlim, add.legend = add.legend, ...)
        if (!is.null(main)) 
            title(main)
    }
    if (type == "I") {
        image.plot(obj$x, obj$y, obj$z, xlab = xlab, ylab = ylab, 
            zlim = zlim, ...)
        if ((!is.null(main)) & type != "b") 
            title(main)
    }
    if (type == "c") {
        if (is.null(levels)) 
            levels <- pretty(obj$z[!is.na(obj$z)], 5)
        contour(obj$x, obj$y, obj$z, levels = levels, labcex = labcex, 
            lwd = 2, ...)
        if ((!is.null(main)) & type != "b") 
            title(main)
    }
    if (type == "b" | type == "C") {
        if (type == "b") {
            par(mar = old.mar)
        }
        image.plot(obj$x, obj$y, obj$z, xlab = xlab, ylab = ylab, 
            graphics.reset = graphics.reset, zlim = zlim, ...)
        if (is.null(levels)) 
            levels <- pretty(obj$z[!is.na(obj$z)], 5)
        contour(obj$x, obj$y, obj$z, add = TRUE, levels = levels, 
            labcex = labcex, col = "black", lwd = 2)
        if ((!is.null(main)) & type != "b") 
            title(main)
    }
    invisible()
}
