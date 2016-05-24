neglogLik <- function(params, object, pmap){
    object <- pmap(object, params)
    return(-logLik(object))
}
