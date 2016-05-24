logLik.gpcm <-
function (object, ...) {
    if (!inherits(object, "gpcm"))
        stop("Use only with 'gpcm' objects.\n")
    out <- object$log.Lik
    df <- sum(sapply(object$coefficients, length))
    if (object$constraint == "1PL")
        df <- df - length(object$coefficients) + 1
    if (object$constraint == "rasch")
        df <- df - length(object$coefficients)    
    attr(out, "df") <- df
    attr(out, "nobs") <- nrow(object$X)
    class(out) <- "logLik"
    out
}
