logLik.sbchoice <- function(object, ...)
{
    out <- summary(object)$loglik
    attr(out, "df") <- sum(!is.na(coefficients(object)))
    attr(out, "nobs") <- object$nobs
    class(out) <- "logLik"
    out
}
