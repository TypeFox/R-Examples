logLik.CNVassoc<-
function (object, ...)
{
    x <- object
    family<-attr(x,"family")
    
    if (family=="binomial") {
        param <- matrix2vector(x$coefficients, x$variant)
        ans <- logLike.logistic(param, x$y, x$X, x$w, x$variant)
    }
    if (family=="poisson") {
        param <- matrix2vector(x$coefficients, x$variant)
        ans <- logLike.poisson(param, x$y, x$X, x$w, x$variant)
    }
    if (family=="gaussian") {
        param <- matrix2vector(x$coefficients, x$variant)
        param <- c(param, x$sigma)
        ans <- logLike.norm(param, x$y, x$X, x$w, x$variant)
    }
    if (family=="weibull") {
        param <- matrix2vector(x$coefficients, x$variant)
        param <- c(param, x$alpha)
        ans <- logLike.weibull(param, x$y, x$cens, x$X, x$w, x$variant)
    }    
    df <- length(param)
    c(logLik = ans, df = df)
}