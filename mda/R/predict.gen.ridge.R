predict.gen.ridge <-
function (object, newdata, ...) 
{
    if (missing(newdata)) 
        fitted(object)
    else scale(newdata, object$xmeans, FALSE) %*% object$coef
}

