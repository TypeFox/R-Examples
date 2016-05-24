print.summary.CNVassoc <-
function (x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    cat("Deviance:", x$deviance, "\n")
    cat("Number of parameters:", x$numparam, "\n")
    cat("Number of individuals:", length(x$y), "\n")    
    cat("\nCoefficients:\n")
    coefs <- x$coefficients
    printCoefmat(coefs, digits = digits, na.print = "", ...)
    family <- attr(x, "family")
    if (family %in% c("binomial","poisson"))
      cat("\n(Dispersion parameter for ", family, " family taken to be ", 1, ")\n\n")
    if (family %in% c("gaussian"))
      cat("\n(Dispersion parameter estimation for ", family, "family is", format(x$sigma^2), ")\n\n")    
    if (family %in% c("weibull"))
      cat("\n(Shape parameter (alpha) estimation for ", family, "family is", format(x$alpha), ")\n\n")    
    if (!family %in% c("gaussian","binomial","poisson","weibull"))
      stop("Family not implemented: it must be 'gaussian', 'binomial', 'poisson' or 'weibull'")
    Varmat <- x$varmat
    p <- NCOL(Varmat)
    cat("\nCovariance between coefficients:\n")
    Varmat <- format(round(Varmat, digits), nsmall = digits,
        digits = digits)
    Varmat[lower.tri(Varmat)] <- ""
    if (family=="gaussian")
        Varmat <- Varmat[-p, -p]
    print(Varmat, quote = FALSE)
    cat("\n")
}
