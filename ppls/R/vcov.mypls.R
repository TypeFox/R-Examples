vcov.mypls=function(object,...){
    dummy<-object$covariance
    if (is.null(dummy)==TRUE){
        cat(paste("WARNING: Covariance of regression coefficients is not available.\n"))
        cat(paste("Returning NULL object.\n"))
    }
    return(dummy)
}
