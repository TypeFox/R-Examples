covmat.default <-
function(x,...) 
{
    object <- list() 
    object$Sigmahat <- x$Sigmahat
    object$lambdahat <- x$lambdahat
    object$Sigmasam <- x$Sigmasam
    object$centered <- x$centered
    object$Target <- x$Target
    class(object) <- "covmat"
    object
}
