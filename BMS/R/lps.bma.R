lps.bma <-
function (object, realized.y, newdata = NULL) 
{
    if (!any(class(object) %in% c("pred.density", "bma", "zlm"))) 
        stop("object must be of class 'pred.density', 'bma' or 'zlm'!")
    if (any(class(object) %in% c("bma", "zlm"))) {
        if (is.null(newdata)) 
            stop("newdata must be provided if object is of class 'bma' or 'zlm'.")
        object = pred.density(object, newdata = newdata)
    }
    return(object$lps(realized.y))
}
