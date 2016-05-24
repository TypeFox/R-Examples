logLik.ltm <-
function (object, ...) {
    if (!inherits(object, "ltm"))
        stop("Use only with 'ltm' objects.\n")
    out <- object$log.Lik
    attr(out, "df") <- if (!is.null(constr <- object$constraint))
        length(object$coef) - nrow(constr) else length(object$coef)
    attr(out, "nobs") <- nrow(object$X)
    class(out) <- "logLik"
    out    
}
