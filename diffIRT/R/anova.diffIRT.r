anova.diffIRT = function (object, object2, ...) {
    if (missing(object2))
        stop("Object 2 is missing.\n")
    if (!(all.equal(object$x, object2$x)))
        warning("It appears that the model in the first object is fitted on a different 
              dataset as compared to the model in the second object.\n")
    LL1 = object$totLL
    LL2 = object2$totLL
    df <- object2$npar - object$npar

    if (df < 0)
        stop("The model in the first object should be nested in the model from the second object.\n")
    LRT <- LL1 - LL2
    if (LRT < 0)
        stop("Problems conducting the likelihood ratio test: Models are not nested, or the solution in object 2 concerns a local minimum.\n")
    p.value = pchisq(LRT, df, lower.tail = FALSE)
    p.value = round(p.value, 3)
    p.value = if (p.value < 0.001) "<0.001" else p.value
    dat <- data.frame(AIC = round(c(object$AIC, object2$AIC), 0), BIC = round(c(object$BIC, object2$BIC), 0),
                sBIC = round(c(object$sBIC, object2$sBIC), 0), DIC = round(c(object$DIC, object2$DIC), 0),
                log.Lik = round(c(LL1, LL2), 3), LRT = c(" ", round(LRT, 2)), df = c("", df),
                p.value = c("", p.value), row.names = c(deparse(substitute(object)), deparse(substitute(object2))))
    cat("\n Likelihood Ratio Table\n")
    print(dat)
    cat("\n")
}
