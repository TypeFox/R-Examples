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
Matern.cor.to.range <- function(d, nu, cor.target = 0.5, 
    guess = NULL, ...) {
    # define local function for root finding
    #
    ftemp <- function(theta, f.extra) {
        Matern(f.extra$d/theta, nu = f.extra$nu) - f.extra$cor.target
    }
    # inital guess is exponential
    if (is.null(guess)) {
        guess[1] <- guess[2] <- -d/log(cor.target)
    }
    #  extra info for function
    f.extra = list(d = d, nu = nu, cor.target = cor.target)
    # find  guesses that are above and below
    while (ftemp(guess[2], f.extra) < 0) {
        guess[2] <- guess[2] * 2
    }
    while (ftemp(guess[1], f.extra) > 0) {
        guess[1] <- guess[1]/2
    }
    temp <- bisection.search(guess[1], guess[2], f = ftemp, f.extra = f.extra, 
        ...)
    return(temp$x)
}
