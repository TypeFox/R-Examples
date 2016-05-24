variable.names.zlm <-
function (object, ...) 
{
    if (!is(object, "zlm")) 
        stop("argument 'object' needs to be zlm object")
    return(names(object$coefficients))
}
