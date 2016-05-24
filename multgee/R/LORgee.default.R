LORgee.default <-
function(x,...) 
{
    object <- list() 
    object$title <- x$title
    object$version <- x$version
    object$link <- x$link
    object$odds.ratio <- x$odds.ratio
    object$terms <- x$terms
    object$contrasts <- x$contrasts
    object$nobs <- x$nobs
    object$convergence <-  x$convergence
    object$coefficients <- x$coefficients
    object$linear.predictors <- x$linear.predictors
    object$fitted.values <- x$fitted.values
    object$residuals <- x$residuals
    object$y <- x$y
    object$id <- x$id
    object$max.id <- x$max.id
    object$clusz <- x$clusz
    object$robust.variance <- x$robust.variance
    object$naive.variance <- x$naive
    object$xnames <- x$xnames
    object$categories <- x$categories
    object$occasions <- x$occasions
    object$LORgee.control <- x$LORgee.control 
    object$ipfp.control <- x$ipfp.control
    object$inverse.method <- x$inverse.method
    object$adding.constant <- x$adding.constant
    object$call <- x$call
    object$trace <- x$trace
    object$pvalue <- x$pvalue
    class(object) <- "LORgee"
    object
}

