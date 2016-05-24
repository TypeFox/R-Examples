
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:                 DESCRIPTION:
#  summary                   Summary method for an object of class 'fREG'
################################################################################


setMethod(f = "summary", signature(object = "fREG"), definition =
    function(object, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Summary method for Regression Modelling, an object of class "fREG"

    # FUNCTION:

    # Digits:
    digits <- max(4, getOption("digits") - 4)

    # Print all from print Method:
    print(object)

    # Add Residual Variance:
    cat("Residual Variance:\n", var(object@fit$residuals))
    cat("\n\n")

    # Internal Function: fResiduals
    fResiduals <-
        function(x, digits)
    {
        cat("Non-Weighted Residuals:\n")
        names = c("Min", "1Q", "Median", "3Q", "Max")
        rq = structure(quantile(x), names = names)
        print(rq, digits = digits)
        names = c("Variance", "StDev", "Skewness", "Kurtosis")
        skewness = sum((x - mean(x))^3/sqrt(var(x))^3)/length(x)
        kurtosis = sum((x - mean(x))^4/var(x)^2)/length(x) - 3
        rq = structure(c(var(x), sqrt(var(x)), skewness, kurtosis),
            names = names)
        print(rq, digits = digits)
        print("done")
        cat("\n")
        invisible()
    }

    # Internal Function: print.summary.LM
    print.summary.LM <- 
        function (x, ...)
    {
        digits = max(4, getOption("digits") - 4)
        symbolic.cor = x$symbolic.cor
        signif.stars = getOption("show.signif.stars")
        # cat("\nCall:\n")
        # cat(paste(deparse(x$call), sep = "\n", collapse = "\n"),
        #   "\n\n", sep = "")
        aliased = x$aliased
        resid = x$residuals
        df = x$df
        rdf = df[2]
        cat(if (!is.null(x$w) && diff(range(x$w)))
            "Weighted ", "Residuals:\n", sep = "")
        if (rdf > 5) {
            nam = c("Min", "1Q", "Median", "3Q", "Max")
            rq = if (length(dim(resid)) == 2)
                structure(apply(t(resid), 1, quantile),
                    dimnames = list(nam, dimnames(resid)[[2]]))
            else structure(quantile(resid), names = nam)
            print(rq, digits = digits, ...)
        } else if (rdf > 0) {
            print(resid, digits = digits, ...)
        } else {
            cat("ALL", df[1], "residuals are 0: no residual ",
                "degrees of freedom!\n")
        }
        if (length(aliased) == 0) {
            cat("\nNo Coefficients\n")
        } else {
            if (nsingular<-df[3] - df[1]) {
                cat("\nCoefficients: (", nsingular, " not defined ",
                    "because of singularities)\n", sep = "")
            } else {
                cat("\nCoefficients:\n")
            }
            coefs = x$coefficients
            if (!is.null(aliased) && any(aliased)) {
                cn = names(aliased)
                coefs = matrix(NA, length(aliased), 4, dimnames =
                    list(cn, colnames(coefs)))
                coefs[!aliased, ] = x$coefficients
            }
            printCoefmat(coefs, digits = digits, signif.stars =
                signif.stars, na.print = "NA", ...)
        }
        cat("\nResidual standard error:", format(signif(x$sigma,
            digits)), "on", rdf, "degrees of freedom\n")
        if (!is.null(x$fstatistic)) {
            cat("Multiple R-Squared:", formatC(x$r.squared,
                digits = digits))
            cat(",  Adjusted R-squared:", formatC(x$adj.r.squared,
                digits = digits), "\nF-statistic:",
                formatC(x$fstatistic[1], digits = digits), "on",
                x$fstatistic[2], "and", x$fstatistic[3],
                "DF,  p-value:", format.pval(pf(x$fstatistic[1],
                x$fstatistic[2], x$fstatistic[3], lower.tail = FALSE),
                digits = digits), "\n") }
        correl = x$correlation
        if (!is.null(correl)) {
            p = NCOL(correl)
            if (p > 1) {
                cat("\nCorrelation of Coefficients:\n")
                if (is.logical(symbolic.cor) && symbolic.cor) {
                    print(symnum(correl, abbr.colnames = NULL)) }
                else {
                    correl = format(round(correl, 2), nsmall = 2,
                      digits = digits)
                    correl[!lower.tri(correl)] = ""
                    print(correl[-1, -p, drop = FALSE], quote = FALSE) }} }
        cat("\n")
        invisible()
    }

    # Internal Function: print.summary.GLM
    print.summary.GLM =
    function (x, ...)
    {
        digits = max(4, getOption("digits") - 4)
        symbolic.cor = x$symbolic.cor
        signif.stars = getOption("show.signif.stars")
        aliased = x$aliased
        df = x$df
        #cat("\nCall:\n")
        #cat(paste(deparse(x$call), sep = "\n", collapse = "\n"),
        #   "\n\n", sep = "")
        cat("Deviance Residuals: \n")
        if (x$df.residual > 5) {
            x$deviance.resid = quantile(x$deviance.resid, na.rm = TRUE)
            names(x$deviance.resid) = c("Min", "1Q", "Median", "3Q",
                "Max") }
        print.default(x$deviance.resid, digits = digits, na.print = "",
            print.gap = 2)
        if (length(aliased) == 0) {
            cat("\nNo Coefficients\n")
        } else {
            if (!is.null(df) && (nsingular = df[3] - df[1]))
                cat("\nCoefficients: (", nsingular, " not defined ",
                "because of singularities)\n", sep = "")
            else cat("\nCoefficients:\n")
            coefs = x$coefficients
            if (!is.null(aliased) && any(aliased)) {
                cn = names(aliased)
                coefs = matrix(NA, length(aliased), 4, dimnames =
                    list(cn, colnames(coefs)))
                coefs[!aliased, ] = x$coefficients }
            printCoefmat(coefs, digits = digits, signif.stars =
                signif.stars, na.print = "NA", ...)
            }
        cat("\n(Dispersion parameter for ", x$family$family,
            " family taken to be ", format(x$dispersion), ")\n\n",
            apply(cbind(paste(format(c("Null", "Residual"),
            width = 8, flag = ""), "deviance:"),
            format(unlist(x[c("null.deviance", "deviance")]),
            digits = max(5, digits + 1)), " on",
            format(unlist(x[c("df.null", "df.residual")])),
            " degrees of freedom\n"), 1, paste, collapse = " "),
            "AIC: ", format(x$aic, digits = max(4, digits + 1)),
            "\n\n", "Number of Fisher Scoring iterations: ",
            x$iter, "\n", sep = "")
        correl = x$correlation
        if (!is.null(correl)) {
            p = NCOL(correl)
            if (p > 1) {
                cat("\nCorrelation of Coefficients:\n")
                if (is.logical(symbolic.cor) && symbolic.cor) {
                    print(symnum(correl, abbr.colnames = NULL)) }
                else {
                    correl = format(round(correl, 2), nsmall = 2,
                      digits = digits)
                    correl[!lower.tri(correl)] = ""
                    print(correl[-1, -p, drop = FALSE], quote = FALSE) }}}
        cat("\n")
        invisible()
    }

    # Internal Function: print.summary.GAM
    print.summary.GAM =
    function(x, ...)
    {
        if (length(x$p.coeff) > 0) {
            cat("Parametric coefficients:\n")
            width = max(nchar(names(x$p.coeff)))
            cat(rep(" ",width), "   Estimate  std. err.    t ratio",
                "    Pr(>|t|)\n", sep = "")
            for (i in 1:length(x$p.coeff))
                cat(formatC(names(x$p.coeff)[i], width = width), " ",
                    formatC(x$p.coeff[i], width=10, digits=5), " ",
                    formatC(x$se[i], width = 10, digits = 4), " ",
                    formatC(x$p.t[i], width = 10, digits = 4), "    ",
                    format.pval(x$p.pv[i]), "\n", sep="") }
        cat("\n")
        if (x$m > 0) {
            cat("Approximate significance of smooth terms:\n")
            width = max(nchar(names(x$chi.sq)))
            cat(rep(" ",width), "        edf       chi.sq     ",
                "p-value\n", sep = "")
            for (i in 1:x$m)
                cat(formatC(names(x$chi.sq)[i], width = width), " ",
                    formatC(x$edf[i], width = 10, digits = 4), "   ",
                    formatC(x$chi.sq[i], width = 10, digits = 5), "     ",
                    format.pval(x$s.pv[i]), "\n", sep = "") }
        cat("\nR-sq.(adj) = ", formatC(x$r.sq, digits = 3, width = 5),
            "   Deviance explained = ", formatC(x$dev.expl*100,
            digits = 3, width = 4), "%", sep = "")
        if (is.null(x$ubre)) {
            cat("\nGCV score = ", formatC(x$gcv, digits = 5), " ", sep = "")
        } else {
            cat("\nUBRE score = ", formatC(x$ubre, digits = 5), sep = "")
        }
        cat("  Scale est. = ", formatC(x$scale, digits = 5,
            width = 8, flag = "-"), "  n = ", x$n, "\n", sep = "")
        invisible()
    }

    # Fit:
    fit <- object@fit

    # Regression Model: LM
    if (object@method == "lm") {
        class(fit) = "lm"
        ans <- stats::summary.lm(object = fit, ...)
        print.summary.LM(x = ans, ...)
    }

    # Regression Model: GLM
    if (object@method == "glm") {
        class(fit) = c("glm", "lm")
        ans <- stats::summary.glm(object = fit, ...)
        print.summary.GLM(x = ans, ...)
    }

    # Regression Model: GAM
    if (object@method == "gam") {
        class(fit) = "gam"
        ans <- mgcv::summary.gam(object = fit, ...)
        print.summary.GAM(x = ans, ...)
    }

    # Regression Model: PPR
    if (object@method == "ppr") {
        # This is what print.ppr produces.
        mu <- fit$mu; ml = fit$ml
        cat("Goodness of fit:\n")
        gof <- fit$gofn; names(gof) = paste(1:ml, "terms")
        print(format(gof[mu:ml], ...), quote = FALSE)
        # This is what summary.ppr produces.
        if (any(fit$edf > 0)) {
            cat("\nEquivalent df for ridge terms:\n")
            edf <- fit$edf
            names(edf) = paste("term", 1:fit$mu)
            print(round(edf, 2), ...)}
    }

    # Regression Model: POLYMARS
    if (object@method == "polymars") {
        class(fit) = "polymars"
        # This is what summary.polymars produces.
        # There is no print.summary.polymars.
        cat("Model Fitting:\n")
        print(fit$fitting)
        if(fit$responses != 1)
            cat("\nResponses:", fit$responses, "\n")
        if(!is.null(fit$Rsquared))
            cat("\nRsquared:",round(fit$Rsquared, 4),"\n")
        cat("\n")
    }

    # Regression Model: NNET
    if (object@method == "nnet") {
        # Use the print Method
    }

    # Return Value:
    invisible()
})


################################################################################
