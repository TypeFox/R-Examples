logLik.zlm <-
function (object, ...) 
{
    if (!is(object, "zlm")) 
        stop("argument 'formula' needs to be zlm object")
    ret = object$marg.lik
    attr(ret, "df") = object$rank + 1
    attr(ret, "nbobs") = object$rank + object$df.residual
    class(ret) = "logLik"
    return(ret)
}
