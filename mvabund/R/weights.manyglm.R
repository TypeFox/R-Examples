weights.manyglm <- function (object, type = c("prior", "working"), ...)
{
    type <- match.arg(type)
    res <- if (type == "prior")
        object$prior.weights
    else object$sqrt.weights
    res <- apply(res, 2, function(x) x^2)
    if (is.null(object$na.action))
        res
    else naresid(object$na.action, res)
}
