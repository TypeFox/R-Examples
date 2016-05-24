predict.polyreg <-
function (object, newdata, ...) 
{
    if (missing(newdata)) {
        z <- fitted(object)
        if (is.null(z)) 
            stop("need to supply newdata")
        else return(z)
    }
    degree <- object$degree
    monomial <- object$monomial
    polybasis(newdata, degree, monomial) %*% object$coef
}

