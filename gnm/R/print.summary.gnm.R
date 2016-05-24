#  Modification of print.summary.glm from the stats package for R.
#
#  Copyright (C) 1995-2006 The R Core Team
#  Copyright (C) 2006, 2008, 2009, 2015 Heather Turner
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

print.summary.gnm <- function (x, digits = max(3, getOption("digits") - 3),
                               signif.stars = getOption("show.signif.stars"),
                               symbolic.cor = x$symbolic.cor, ...)
{
    cat("\nCall:\n", deparse(x$call), "\n", sep = "", fill = TRUE)

    cat("Deviance Residuals: \n")
    if (length(x$deviance.resid) > 5) {
        x$deviance.resid <- quantile(x$deviance.resid, na.rm = TRUE)
        names(x$deviance.resid) <- c("Min", "1Q", "Median", "3Q",
            "Max")
    }
    print.default(x$deviance.resid, digits = digits, na.print = "",
                  print.gap = 2)

    tidy.zeros <- function(vec)
        ifelse(abs(vec) < 100 * .Machine$double.eps, 0, vec)
    coefs <- tidy.zeros(coef(x))
    if (length(ofInterest(x)) > 0)
        coefs <- coefs[ofInterest(x), , drop = FALSE]
    non.elim <- length(coefs)
    elim <- length(x$eliminated)

    if (non.elim | elim) {
        cat("\nCoefficients", " of interest"[!is.null(ofInterest(x))], ":\n",
            sep = "")
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
            signif.legend = !elim, na.print = "NA", ...)
        if (elim){
            cat("\nEliminated coefficients:\n", sep = "")
            printCoefmat(x$eliminated, digits = digits,
                         signif.stars = signif.stars, na.print = "NA", ...)
        }
        coefs <- c(coefs[,2], x$eliminated[,2])
        if (any(!is.na(coefs)))
            cat("\n(Dispersion parameter for ", x$family$family,
                " family taken to be ", format(x$dispersion), ")\n", sep = "")
        if (any(is.na(coefs)))
            cat("\nStd. Error is NA where coefficient has been constrained or",
                "is unidentified\n")
    }
    else cat("\nNo coefficients", " of interest"[!is.null(ofInterest(x))],
             ". \n\n", sep = "")

    cat("\nResidual deviance: ", format(x$deviance,
                                        digits = max(5, digits + 1)),
        " on ", format(x$df.residual, digits = max(5, digits + 1)),
        " degrees of freedom\n",
        "AIC: ", format(x$aic, digits = max(4, digits + 1)), "\n\n",
        "Number of iterations: ", x$iter, "\n", sep = "")
    correl <- x$correlation
    if (!is.null(correl)) {
        if (attr(x$cov.scaled, "eliminate")) {
            eliminate <- seq(attr(x$cov.scaled, "eliminate"))
            correl <- correl[-eliminate, -eliminate]
        }
        p <- NCOL(correl)
        if (p > 1) {
            cat("\nCorrelation of Coefficients:\n")
            if (is.logical(symbolic.cor) && symbolic.cor) {
                print(symnum(correl, abbr.colnames = NULL))
            }
            else {
                correl <- format(round(correl, 2), nsmall = 2, digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop = FALSE], quote = FALSE)
            }
        }
    }
    cat("\n")
    invisible(x)
}
