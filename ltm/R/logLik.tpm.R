logLik.tpm <-
function (object, ...) {
    if (!inherits(object, "tpm"))
        stop("Use only with 'tpm' objects.\n")
    out <- object$log.Lik
    attr(out, "df") <- if (!is.null(constr <- object$constraint)) {
        if (object$type == "rasch") 2 * nrow(object$coef) + 1 - nrow(constr) else 3 * nrow(object$coef) - nrow(constr)
    } else {
        if (object$type == "rasch") 2 * nrow(object$coef) + 1 else 3 * nrow(object$coef)
    }
    attr(out, "nobs") <- nrow(object$X)
    class(out) <- "logLik"
    out
}
