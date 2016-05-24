logLik.dbchoice <- function(object, ...)
{
    out <- object$loglik
    attr(out, "df") <- sum(!is.na(coefficients(object)))
    attr(out, "nobs") <- object$nobs
    class(out) <- "logLik"
    out
}
