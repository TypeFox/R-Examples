##################################################################
## L2 Group family methods
##################################################################
setMethod("LogDeriv", signature(object = "L2GroupParamFamily"),
           function(object) object@LogDeriv)

setMethod("locscalename", signature(object = "L2LocationScaleUnion"),
           function(object) object@locscalename)

setMethod("scaleshapename", signature(object = "L2ScaleShapeUnion"),
           function(object) object@scaleshapename)

setMethod("scalename", signature(object = "L2LocationScaleUnion"),
           function(object) object@locscalename["scale"])

setMethod("scalename", signature(object = "L2ScaleShapeUnion"),
           function(object) object@scaleshapename["scale"])

setMethod("withPosRestr", signature(object = "L2ScaleShapeUnion"),
           function(object) object@param@withPosRestr)

setReplaceMethod("LogDeriv", "L2GroupParamFamily",
    function(object, value){
        object@LogDeriv <- value
        object
    })

setReplaceMethod("locscalename", "L2LocationScaleUnion",
    function(object, value){
        if(length(value)>2||length(value)<1)
           stop("value of slot 'locscalename' must be of length one or two")
        if(is.null(names(value))) names(value) <- c("loc","scale")
        object@locscalename <- value
        object
    })

setReplaceMethod("scaleshapename", "L2ScaleShapeUnion",
    function(object, value){
        if(length(value)>2||length(value)<1)
           stop("value of slot 'scaleshapename' must be of length one or two")
        if(is.null(names(value))) names(value) <- c("scale","shape")
        object@scaleshape <- value
        object
    })

setReplaceMethod("withPosRestr", "L2ScaleShapeUnion",
    function(object, value){
        if(length(value)!=1)
           stop("value of slot 'withPosRestr' must be of length one")
        param <- object@param
        withPosRestr(param) <- value
        object@param <- param
        object
    })

