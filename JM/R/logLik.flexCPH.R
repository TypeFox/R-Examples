logLik.flexCPH <-
function (object, ...) {
    out <- object$logLik
    attr(out, "df") <- length(unlist(object$coefficients))
    attr(out, "nobs") <- length(object$d)
    class(out) <- "logLik"
    out    
}
