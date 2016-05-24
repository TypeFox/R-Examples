setGeneric("nchannel",
function(object) standardGeneric("nchannel"))

setMethod("nchannel", signature(object = "Wave"), 
function(object){
    object@stereo + 1
})

setMethod("nchannel", signature(object = "WaveMC"), 
function(object){
    ncol(object)
})
