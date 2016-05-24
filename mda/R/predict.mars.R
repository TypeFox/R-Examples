predict.mars <-
function (object, newdata, ...) 
{
    if (missing(newdata)) {
        z <- fitted(object)
        if (is.null(z)) 
            stop("need to supply newdata")
        else return(z)
    }
    model.matrix.mars(object, newdata) %*% object$coefficients
}

