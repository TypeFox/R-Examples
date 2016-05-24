summary.mpp <- function(object, ...){
    object$LL <- logLik(object)
    object$data <- NULL
    if (inherits(object$gif, "function")) object$gif <- NULL
    class(object) <- c("summary.mpp", "list")
    return(object)
}

