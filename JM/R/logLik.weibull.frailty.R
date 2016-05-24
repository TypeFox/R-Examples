logLik.weibull.frailty <-
function (object, ...) {
    out <- object$logLik
    attr(out, "df") <- length(unlist(object$coefficients))
    attr(out, "nobs") <- length(unique(object$id))
    class(out) <- "logLik"
    out    
}
