#  Copyright (C) 2005 Heather Turner
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

psum <- function(...) {
    summandList <- list(...)
    nSummands <- length(summandList)
    if (nSummands == 0) return(0)
    else if (nSummands == 1) return(summandList[[1]])
    else {
        trySum <- try(summandList[[1]] + do.call("Recall", summandList[-1]),
                      silent = TRUE)
        if (inherits(trySum, "try-error"))
            stop("addition not implemented for types of argument supplied")
        else trySum
    }
}
