deviance.bma <-
function (object, exact = FALSE, ...) 
{
    if (is.bma(object)) {
        xx = as.matrix(object$X.data)
        ebeta = estimates.bma(object, order.by.pip = FALSE, exact = exact)[, 
            2, drop = TRUE]
    }
    else if (is(object, "lm")) {
        xx = as.matrix(object$model)
        ebeta = coef(object)
        if (length(ebeta) == ncol(xx)) 
            ebeta = ebeta[-1]
    }
    else stop("Required input is an object of class 'bma' or 'lm'/'zlm'.")
    xx = xx - matrix(colMeans(xx), nrow(xx), ncol(xx), byrow = TRUE)
    ess = as.vector(crossprod(ebeta, as.vector(crossprod(xx[, 
        -1, drop = FALSE], xx[, 1]))))
    return((as.vector(crossprod(xx[, 1, drop = TRUE])) - ess))
}
