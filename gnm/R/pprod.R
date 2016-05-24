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

pprod <- function(...) {
    factorList <- list(...)
    nFactors <- length(factorList)
    if (nFactors == 0) return(1)
    else if (nFactors == 1) return(factorList[[1]])
    else {
        tryProduct <- try(factorList[[1]] * do.call("Recall", factorList[-1]),
                          silent = TRUE)
        if (inherits(tryProduct, "try-error"))
            stop("multiplication not implemented for types of argument supplied")
        else tryProduct
    }
}

