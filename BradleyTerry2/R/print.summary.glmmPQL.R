print.summary.BTglmmPQL <- function(x, digits = max(3, getOption("digits") - 3),
                                   symbolic.cor = x$symbolic.cor,
                                   signif.stars = getOption("show.signif.stars"),
                                   ...)  {
    cat("\nCall:\n", deparse(x$call), sep = "", fill = TRUE)
    p <- length(x$aliased)
    tidy.zeros <- function(vec) ifelse(abs(vec) < 100 * .Machine$double.eps,
                                       0, vec)
    if (p == 0) {
        cat("\nNo Fixed Effects\n")
    }
    else {
        if (nsingular <- p - x$rank) {
            cat("\nFixed Effects: (", nsingular,
                " not defined because of singularities)\n", sep = "")
            cn <- names(x$aliased)
            pars <- matrix(NA, p, 4, dimnames = list(cn, colnames(x$fixef)))
            pars[!x$aliased, ] <- tidy.zeros(x$fixef)
        }
        else {
            cat("\nFixed Effects:\n")
            pars <- tidy.zeros(x$fixef)
        }
        printCoefmat(pars, digits = digits, signif.stars = signif.stars,
            na.print = "NA", ...)
    }
    cat("\n(Dispersion parameter for ", x$family$family,
        " family taken to be 1)\n", sep = "")
    cat("\nRandom Effects:\n")
    pars <- tidy.zeros(x$ranef)
    printCoefmat(pars, digits = digits, signif.stars = signif.stars,
            na.print = "NA", ...)
    if (nzchar(mess <- naprint(x$na.action)))
        cat("\n", mess, "\n", sep = "")
    cat("\nNumber of iterations: ", x$iter, "\n", sep = "")
    correl <- x$correlation
    if (!is.null(correl)) {
        if (x$rank > 1) {
            cat("\nCorrelation of Coefficients:\n")
            if (is.logical(symbolic.cor) && symbolic.cor) {
                print(symnum(correl, abbr.colnames = NULL))
            }
            else {
                correl <- format(round(correl, 2), nsmall = 2,
                  digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -x$rank, drop = FALSE], quote = FALSE)
            }
        }
    }
    cat("\n")
    invisible(x)
}
