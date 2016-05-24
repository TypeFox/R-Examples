deviance.rma <-
function (object, REML, ...) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    if (missing(REML)) {
        if (object$method == "REML") {
            REML <- TRUE
        }
        else {
            REML <- FALSE
        }
    }
    if (REML) {
        val <- object$fit.stats$REML[2]
    }
    else {
        val <- object$fit.stats$ML[2]
    }
    return(val)
}
