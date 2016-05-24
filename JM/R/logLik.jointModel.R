logLik.jointModel <-
function (object, ...) {
    if (!inherits(object, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    out <- object$logLik
    attr(out, "df") <- nrow(object$Hessian)
    attr(out, "nobs") <- object$n
    class(out) <- "logLik"
    out
}
