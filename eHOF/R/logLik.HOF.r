"logLik.HOF" <-
    function(object, ...)
{
     sapply(object$models, function(x) x$deviance/(-2))
}
