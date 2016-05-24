setGeneric("equalWave",
function(object1, object2) standardGeneric("equalWave"))

setMethod("equalWave", signature(object1 = "Wave", object2 = "Wave"), 
function(object1, object2){
    if(!(validObject(object1) && validObject(object2)))
        stop("Not a valid 'Wave' object")
    if(object1@samp.rate != object2@samp.rate)
        stop("Sampling Rate of the 'Wave' objects differ")    
    if(object1@bit != object2@bit)
        stop("Bit resolution of the 'Wave' objects differ")    
    if(object1@pcm != object2@pcm)
        stop("Format (pcm or float) of the 'Wave' objects differ")    
    if(xor(object1@stereo, object2@stereo))
        stop("One 'Wave' object is mono, the other one stereo")
}
)


setMethod("equalWave", signature(object1 = "WaveMC", object2 = "WaveMC"), 
function(object1, object2){
    if(!(validObject(object1) && validObject(object2)))
        stop("Not a valid 'WaveMC' object")
    if(object1@samp.rate != object2@samp.rate)
        stop("Sampling Rate of the 'WaveMC' objects differ")    
    if(object1@bit != object2@bit)
        stop("Bit resolution of the 'WaveMC' objects differ")    
    if(object1@pcm != object2@pcm)
        stop("Format (pcm or float) of the 'WaveMC' objects differ")    
    if(ncol(object1) != ncol(object2))
        stop("The number of channels of the 'WaveMC' objects differ")
}
)


setMethod("equalWave", signature(object1 = "WaveMC", object2 = "ANY"), 
function(object1, object2){
    stop("Both arguments have to be either of class 'Wave' or 'WaveMC'.")
}
)


setMethod("equalWave", signature(object1 = "Wave", object2 = "ANY"), 
function(object1, object2){
    stop("Both arguments have to be either of class 'Wave' or 'WaveMC'.")
}
)
