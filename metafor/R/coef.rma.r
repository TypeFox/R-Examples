coef.rma <-
function (object, ...) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    x <- object
    coefs <- c(x$b)
    names(coefs) <- rownames(x$b)
    return(coefs)
}
