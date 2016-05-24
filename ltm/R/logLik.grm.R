logLik.grm <-
function (object, ...) {
    if (!inherits(object, "grm"))
        stop("Use only with 'grm' objects.\n")
    out <- object$log.Lik
    df <- sapply(object$coef, length)
    attr(out, "df") <- if (object$constrained) sum(df) - length(df) + 1 else sum(df)
    attr(out, "nobs") <- nrow(object$X)
    class(out) <- "logLik"
    out
}
