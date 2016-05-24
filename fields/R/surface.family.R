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
"surface.Krig" <- function(object, grid.list = NULL, 
    extrap = FALSE, graphics.reset = NULL, xlab = NULL, ylab = NULL, 
    main = NULL, zlab = NULL, zlim = NULL, levels = NULL, type = "C", 
    nx = 80, ny = 80, ...) {
    ## modified so that you can give main, and ylab as arguments
    ## in ... and have them passed correctly
    out.p <- predictSurface(object, grid.list = grid.list, extrap = extrap, 
        nx = nx, ny = ny, drop.Z = TRUE)
    if (!is.null(ylab)) 
        out.p$ylab <- ylab
    if (!is.null(xlab)) 
        out.p$xlab <- xlab
    if (!is.null(zlab)) 
        out.p$zlab <- zlab
    if (!is.null(main)) 
        out.p$main <- main
    ##    else
    ##      out.p$main <- NULL
    plot.surface(out.p, type = type, graphics.reset = graphics.reset, 
        levels = levels, zlim = zlim, ...)
    invisible()
}
# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"surface" <- function(object, ...) {
    UseMethod("surface")
}
# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"surface.default" <- function(object, ...) {
    plot.surface(object, ...)
}
# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"surface.mKrig" <- function(object, grid.list = NULL, 
    extrap = FALSE, graphics.reset = NULL, xlab = NULL, ylab = NULL, 
    main = NULL, zlab = NULL, zlim = NULL, levels = NULL, type = "C", 
    nx = 80, ny = 80, ...) {
    ## modified so that you can give main, and ylab as arguments
    ## in ... and have them passed correctly
    out.p <- predictSurface(object, grid.list = grid.list, extrap = extrap, 
        nx = nx, ny = ny, drop.Z = TRUE)
    if (!is.null(ylab)) 
        out.p$ylab <- ylab
    if (!is.null(xlab)) 
        out.p$xlab <- xlab
    if (!is.null(zlab)) 
        out.p$zlab <- zlab
    if (!is.null(main)) 
        out.p$main <- main
    ##    else
    ##      out.p$main <- NULL
    plot.surface(out.p, type = type, graphics.reset = graphics.reset, 
        levels = levels, zlim = zlim, ...)
    invisible()
}
# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#"surface.surface" <- function(object, ...) {
#    #
#    plot.surface(object, ...)
#}
