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
"bisection.search" <- function(x1, x2, f, tol = 1e-07, 
    niter = 25, f.extra = NA, upcross.level = 0) {
    f1 <- f(x1, f.extra) - upcross.level
    f2 <- f(x2, f.extra) - upcross.level
    if (f1 > f2) 
        stop(" f1 must be < f2 ")
    iter <- niter
    for (k in 1:niter) {
        xm <- (x1 + x2)/2
        fm <- f(xm, f.extra) - upcross.level
        if (fm < 0) {
            x1 <- xm
            f1 <- fm
        }
        else {
            x2 <- xm
            f2 <- fm
        }
        if (abs(fm) < tol) {
            iter <- k
            break
        }
    }
    xm <- (x1 + x2)/2
    fm <- f(xm, f.extra) - upcross.level
    list(x = xm, fm = fm, iter = iter)
}
