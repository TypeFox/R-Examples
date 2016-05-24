#  Copyright (C) 2010-2012 Heather Turner
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

print.meanResiduals <- function (x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nModel call:\n", deparse(attr(x, "call"),
                                   width.cutoff = options()$width),
        "\n", sep = "", fill = TRUE)

    cat("Mean residuals by ", attr(x, "by"),  ":\n\n", sep = "")
    if (!inherits(x, "table")) x <- as.numeric(x)
    NextMethod(object = x, digits = digits, print.gap = 2, ...)
}
