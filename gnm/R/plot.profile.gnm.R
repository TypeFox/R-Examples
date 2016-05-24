#  Modification of plot.profile from the stats package for R.
#
#  File MASS/profiles.q copyright (C) 1996 D. M. Bates and W. N. Venables.
#  port to R by B. D. Ripley copyright (C) 1998
#  corrections copyright (C) 2000,3,6 B. D. Ripley
#  Copyright (C) 2006 Heather Turner
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

plot.profile.gnm <- function (x, nseg, ...)
{
    nulls <- sapply(x, is.null)
    if (all(nulls))
        return(NULL)
    x <- x[!nulls]
    pnames <- names(x)
    pnames <- pnames[!is.na(x[pnames])]
    nr <- ceiling(sqrt(length(pnames)))
    oldpar <- par(mfrow = c(nr, nr))
    on.exit(par(oldpar))
    for (nm in pnames) {
        z <- x[[nm]][[1]]
        parval <- x[[nm]][[2]][, nm]
        plot(parval, z, xlab = nm, ylab = "z", type = "n")
        if (sum(z == 0) == 1)
            points(parval[z == 0], 0, pch = 3)
        splineVals <- spline(parval, z)
        lines(splineVals$x, splineVals$y)
    }
}
