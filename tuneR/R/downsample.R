downsample <- 
function(object, samp.rate){
    if(!(is(object, "Wave") || is(object, "WaveMC"))) 
        stop("'object' needs to be of class 'Wave' or 'WaveMC'")
    validObject(object)
    if((!is.numeric(samp.rate)) || (samp.rate < 2000) || (samp.rate > 192000))
            stop("samp.rate must be an integer in [2000, 192000].")
    if(object@samp.rate > samp.rate){
        ll <- length(object)
        object <- object[seq(1, ll, length = samp.rate * ll / object@samp.rate)]
        object@samp.rate <- samp.rate
    }
    else warning("samp.rate > object's original sampling rate, hence object is returned unchanged.")
    return(object)
}
