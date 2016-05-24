"coef.drc" <-
function(object, ...)
{
    if (!is.null(object$"coefficients"))
    {
        return(object$"coefficients")
    } else {
        retVec <- object$fit$par
        names(retVec) <- object$parNames[[1]]
        return(retVec)
    }
}
