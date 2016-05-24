#  Copyright (C) 2005-2008, 2010 Heather Turner and David Firth
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

print.gnm <- function (x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n", deparse(x$call), "\n", sep = "", fill = TRUE)
    if (length(coef(x)) && (is.null(ofInterest(x)) || length(ofInterest(x)))) {
        cat("Coefficients", " of interest"[!is.null(ofInterest(x))], ":\n",
            sep = "")
        if (!is.null(ofInterest(x)))
            print.default(format(coef(x)[ofInterest(x)], digits = digits),
                          print.gap = 2, quote = FALSE)
        else
            print.default(format(coef(x), digits = digits), print.gap = 2,
                          quote = FALSE)
    }
    else cat("No coefficients", " of interest"[!is.null(ofInterest(x))],
             ". \n\n", sep = "")
    cat("\nDeviance:           ", format(x$deviance, digits),
        "\nPearson chi-squared:",
        format(sum(na.omit(c(residuals(x, type = "pearson")))^2), digits),
        "\nResidual df:        ", x$df.residual, "\n")
    invisible(x)
}
