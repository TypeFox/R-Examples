variable.names.bma <-
function (object, ...) 
{
    if (!is.bma(object)) 
        stop("argument 'object' needs to be a bma object")
    return(c("(Intercept)", object$reg.names))
}
