logLik.rasch <-
function (object, ...) {
    if (!inherits(object, "rasch"))
        stop("Use only with 'rasch' objects.\n")
    out <- object$log.Lik
    attr(out, "df") <- if (!is.null(constr <- object$constraint))
        nrow(object$coef) + 1 - nrow(constr) else nrow(object$coef) + 1
    attr(out, "nobs") <- nrow(object$X)
    class(out) <- "logLik"
    out    
}
