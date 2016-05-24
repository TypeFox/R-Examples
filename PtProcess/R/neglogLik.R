neglogLik <- function(params, object, pmap=NULL, SNOWcluster=NULL){
    #   have an arg to subset parameters
    if (is.null(pmap)) object$params <- params
    else object <- pmap(object, params)
    x <- -logLik(object, SNOWcluster)
    if (is.infinite(x) | is.na(x)) return(1e15)
    else return(x)
}

