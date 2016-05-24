setMethod("length", signature(x = "Wave"), 
function(x){
    validObject(x)
    length(x@left)
})

setMethod("length", signature(x = "WaveMC"), 
function(x){
    validObject(x)
    nrow(x@.Data)
})
