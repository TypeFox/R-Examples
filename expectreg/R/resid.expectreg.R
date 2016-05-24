resid.expectreg <-
function (object, ...) 
{
    fit.res = object$response - object$fitted
    val.res = list()
    for (i in 1:length(object$values)) val.res[[i]] = object$response - 
        object$values[[i]]
    list(fitted = fit.res, values = val.res)
}
