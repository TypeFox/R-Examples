logLik.survreg <-
function (object, ...)
{
    val <- object$loglik
    attr(val, "df") <- object$df
    class(val) <- "logLik"
    val
}
