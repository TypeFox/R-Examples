setGeneric("bind",
function(object, ...) standardGeneric("bind"))

setMethod("bind", signature(object = "Wave"), 
function(object, ...){
    allobjects <- as.list(list(...))
    lapply(allobjects, equalWave, object)
    allobjects <- c(list(object), allobjects)
    object@left <- unlist(lapply(allobjects, slot, "left"))
    if(object@stereo)
        object@right <- unlist(lapply(allobjects, slot, "right"))
    return(object)
}
)

setMethod("bind", signature(object = "WaveMC"), 
function(object, ...){
    allobjects <- as.list(list(...))
    lapply(allobjects, equalWave, object)
    allobjects <- c(list(object), allobjects)
    object@.Data <- do.call(rbind, allobjects)
    return(object)
}
)
