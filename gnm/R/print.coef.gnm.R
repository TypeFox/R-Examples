#  Copyright (C) 2005, 2006, 2010 Heather Turner
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

print.coef.gnm <- function(x, ...) {
    if (!is.null(attr(x, "ofInterest"))) {
        if (length(attr(x, "ofInterest"))){
            cat("Coefficients of interest:\n", sep = "")
            print.default(format(x[attr(x, "ofInterest")]), quote = FALSE)
        }
        else
            cat("No coefficients of interest\n")
    }
    else {
        cat("Coefficients:\n")
        print.default(format(x), quote = FALSE)
    }
}
