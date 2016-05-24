#  as tolNorm2 method in rankMatrix from the Matrix package, but avoids validity
#  checks - much faster if need to do repeated rank calculations
#
#  Copyright (C) 2007 Martin Maechler
#  Copyright (C) 2010 Heather Turner
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

quickRank <- function(X, tol = NULL) {
    sval <- svd(X, 0, 0)$d
    if (is.null(tol))
        sum(sval >= max(dim(X)) * .Machine$double.eps * sval[1])
    else
        sum(sval >= tol)
}
