predict.TML <- function(object, newdata = NULL, ...) 
{
    beta <- object$th1
    if (missing(newdata) || is.null(newdata))
        pred <- fitted.values(object)
    else {
        if (is.vector(newdata))
            X <- newdata
        else X <- as.matrix(newdata)
        pred <- X %*% beta
    }
    pred
}


    
    