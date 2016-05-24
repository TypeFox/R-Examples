
predict.coxr <- function(object, newdata, ...) {

    if ( !inherits(object, "coxr") ) {
        stop("use only with \"coxr\" objects")
    }
    
    beta = object$coef
    if ( length(object$skip) > 0 ) {
        beta[object$skip] = 0
    }

    if ( missing(newdata) ) {
        predictor <- as.vector(object$x %*% beta)
    } else {
        predictor <- as.vector(newdata %*% beta)
    }

    predictor

}