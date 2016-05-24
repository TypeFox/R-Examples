chisq.irt <-
function (object) {
    if (!class(object) %in% c("ltm", "rasch", "tpm"))
        stop("Use only with 'ltm', 'rasch' or 'tpm' objects.\n")
    nam.obj <- deparse(substitute(object))
    p <- ncol(object$X)
    X <- expand.grid(rep(list(c(0, 1)), p))
    res <- residuals(object, resp.patterns = data.matrix(X))
    Tstat <- sum((res[, "Resid"])^2)
    df <- 2^p - attr(logLik(object), "df") - 1
    pval <- pchisq(Tstat, df, lower.tail = FALSE)
    names(Tstat) <- "X^2"
    names(df) <- "df"
    rval <- list(statistic = Tstat, parameter = df, p.value = pval,
        alternative = "the fitted IRT model is not appropriate for the data set",
        method = "Pearson X2 for IRT Models",
        data.name = paste("data from model '", nam.obj, "'", sep = ""))
    class(rval) <- "htest"
    rval
}
