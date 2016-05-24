"fitted.drc" <-
function(object, ...)
{
#    if (missing(...))
#    {
#        return(object$"predres"[, 1])
#    } else {
#        predict(object, ...)
#    }
    predict(object, ...)
##    return(object$"predres"[, 1])
}
