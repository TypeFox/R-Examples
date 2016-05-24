logLik.tegarch <-
function(object, ...)
{
out <- object$objective
attr(out, "nobs") <- length(object$y)
attr(out, "df") <- length(object$par)
class(out) <- "logLik"
return(out)
}
