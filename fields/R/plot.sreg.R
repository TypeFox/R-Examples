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
"plot.sreg" <- function(x, digits = 4, which = 1:4, 
    ...) {
    out <- x
    if (any(which == 1)) {
        plot(out$x, out$y, ylab = "predicted", xlab = " X", bty = "n", 
            ...)
        matlines(out$predicted$x, out$predicted$y, lty = 1)
    }
    if (any(which == 2) & length(out$lambda) == 1) {
        plot(out$fitted.values, out$residuals, ylab = "residuals", 
            xlab = " predicted values", bty = "n", ...)
        yline(0)
    }
    if (any(which == 3)) {
        if (nrow(out$gcv.grid) > 1) {
            # trim off + infinity due to pole in the denominator of GCV function
            #with cost
            ind <- out$gcv.grid[, 3] < 1e+19
            out$gcv.grid <- out$gcv.grid[ind, ]
            yr <- range(unlist(out$gcv.grid[, 3:5]), na.rm = TRUE)
            plot(out$gcv.grid[, 2], out$gcv.grid[, 3], xlab = "Eff. parameters", 
                ylab = " GCV function", bty = "n", ylim = yr, 
                log = "y", ...)
            lines(out$gcv.grid[, 2], out$gcv.grid[, 4], lty = 2)
            lines(out$gcv.grid[, 2], out$gcv.grid[, 5], lty = 1)
            xline(out$eff.df)
            title("GCV-points , solid- GCV model,\ndashed- GCV one", 
                cex = 0.6)
        }
    }
    if (any(which == 4)) {
        if (length(out$lambda) == 1) {
            hist(out$residuals, xlab = "Residuals", main = "")
        }
        else {
            bplot(out$residuals, names = format(round(out$trace, 
                1)), xlab = "eff df")
            title("Residuals")
        }
    }
}
