#' @rdname pcoa.sig
#' @encoding UTF-8
#' @export
summary.pcoasig<-function(object, ...){
    res<-list()
    res$call<-object$call
    res$values<-object$PCoA$values
    res$vectors<-object$PCoA$vectors
    res$correlations<-object$correlations
    res$probabilities<-object$probabilities
    return(res)
}