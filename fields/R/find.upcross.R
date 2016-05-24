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
"find.upcross" <- function(fun, fun.info, upcross.level = 0, 
    guess = 1, tol = 1e-05) {
    l1 <- guess
    tr <- 0
    for (k in 1:50) {
        tr <- fun(l1, fun.info) - upcross.level
        if (tr >= 0) 
            break
        else {
            guess <- l1
        }
        l1 <- l1 * 2
    }
    if (tr < 0) {
        warning("Failed to find the upcrossing")
        return(NA)
    }
    tr <- 0
    l2 <- guess
    for (k in 1:50) {
        tr <- fun(l2, fun.info) - upcross.level
        if (tr <= 0) 
            break
        l2 <- l2/2
    }
    if (tr > 0) {
        warning("Failed to find the upcrossing")
        return(NA)
    }
    out <- bisection.search(l2, l1, fun, tol = tol, f.extra = fun.info, 
        upcross.level = upcross.level)$x
    (out)
}
