#' @rdname pcps
#' @encoding UTF-8
#' @export
summary.pcps<-function(object, ...){
    res<-list()
    res$values<-object$values
    res$vectors<-object$vectors
    res$correlations<-object$correlations
    return(res)
}