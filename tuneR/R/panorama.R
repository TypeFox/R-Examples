setGeneric("panorama",
function(object, pan = 1) standardGeneric("panorama"))

setMethod("panorama", signature(object = "Wave"), 
function(object, pan = 1){
    validObject(object)
    if(!is.numeric(pan) || abs(pan) > 1)
        stop("'pan' must be numeric in [-1, 1].")
    if(!object@stereo){
        warning("panorama called on a mono Wave object, returned object is unchanged")
        return(object)
    }
    pan <- pan/2
    right <- object@right
    object@right <- (0.5 + pan) * right + (0.5 - pan) * object@left
    object@left  <- (0.5 - pan) * right + (0.5 + pan) * object@left
    return(object)
}
)

setMethod("panorama", signature(object = "WaveMC"), 
function(object, pan = 1){
    if(ncol(object) > 2)
        stop("object needs to be a stereo Wave object")
    validObject(object)
    if(!is.numeric(pan) || abs(pan) > 1)
        stop("'pan' must be numeric in [-1, 1].")
    if(ncol(object)==1){
        warning("panorama called on a mono WaveMC object, returned object is unchanged")
        return(object)
    }
    pan <- pan/2
    right <- object@.Data[,2]
    object@.Data[,2] <- (0.5 + pan) * right + (0.5 - pan) * object@.Data[,1]
    object@.Data[,1]  <- (0.5 - pan) * right + (0.5 + pan) * object@.Data[,1]
    return(object)
}
)
