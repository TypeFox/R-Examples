anova.ltm <-
function (object, object2, ...) {
    if (!inherits(object, "ltm"))
        stop("Use only with 'ltm' objects.\n")
    if (missing(object2))
        stop("anova.ltm() computes LRTs between two nested IRT models.")
    if (!inherits(object2, "ltm") && !inherits(object2, "tpm"))
        stop("'object2' must inherit from either class 'ltm' or class 'tpm'.\n")
    if (!isTRUE(all.equal(object$X, object2$X)))
        warning("it seems that the two models are fitted in different data sets.")
    L0 <- logLik(object)
    L1 <- logLik(object2)
    nb0 <- attr(L0, "df")
    nb1 <- attr(L1, "df")
    df. <- nb1 - nb0
    if (df. < 0)
        stop("'object' is not nested in 'object2'.\n")
    LRT <- - 2 * (L0 - L1)
    attributes(LRT) <- NULL
    if (LRT < 0)
        warning("either the two models are not nested or the model represented by 'object2' fell on a local maxima.\n")
    p.value <- pchisq(LRT, df., lower.tail = FALSE)
    out <- list(nam0 = deparse(substitute(object)), L0 = L0, aic0 = AIC(object), 
                bic0 = AIC(object, k = log(attr(L0, "n"))), nam1 = deparse(substitute(object2)), L1 = L1, 
                aic1 = AIC(object2), bic1 = AIC(object2, k = log(attr(L1, "n"))), LRT = LRT, df = df., 
                p.value = p.value, call = object$call)
    class(out) <- "aov.ltm"
    out
}
