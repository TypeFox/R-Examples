#  Copyright (C) 2006, 2008, 2009 Heather Turner
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

print.profile.gnm <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
    #if (attr(x, "eliminate"))
       # coefs <- coefs[-seq(attr(x$cov.scaled, "eliminate")), ]

    if (length(x)) {
        if (any(sapply(x, function(x) isTRUE(is.na(x)))))
            cat("\nProfile is NA where coefficient has been constrained or",
                "is unidentified\n\n")
        print.default(x)
    }
    else cat("\nNo coefficients profiled.\n\n", sep = "")
    invisible(x)
}
