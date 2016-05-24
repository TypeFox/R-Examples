##########
# define class Wave
setClass("WaveGeneral",
    slots = representation(samp.rate = "numeric", bit = "numeric", pcm = "logical"),
    prototype = prototype(samp.rate = 44100, bit = 16, pcm = TRUE))

setClass("Wave",
    slots = representation(left = "numeric",
    right = "numeric", stereo = "logical"),
    contains = "WaveGeneral",
    prototype = prototype(stereo = TRUE))

setClass("WaveMC", contains=c("WaveGeneral", "matrix"))


## convert Wave objects from tuneR <= 0.4-1 to the extended representation
updateWave <- function(object){
     if(!.hasSlot(object, "pcm"))
        object@pcm <- TRUE
     object
}

setValidity("WaveGeneral", 
function(object){
    if(!((length(object@samp.rate) < 2) && (object@samp.rate > 0)))
            return("slot 'samp.rate' of a Wave object must be a positive numeric of length 1")
    if(!((length(object@bit) < 2) && (object@bit %in% c(1, 8, 16, 24, 32, 64))))
            return("slot 'bit' of a Wave object must be a positive numeric (1, 8, 16, 24, 32 or 64) of length 1")
    if(!((length(object@pcm) < 2)))
        return("slot 'pcm' of a Wave object must be a logical of length 1")
    if(object@pcm && object@bit==64)
        return("pcm Wave objects must have a resolution < 64 bit")
    if((!object@pcm) && !(object@bit %in% c(32, 64)))
        return("Float (non pcm) Wave objects must have a resolution of either 32 or 64 bit")
    return(TRUE)
})

setValidity("Wave", 
function(object){
    if(!(length(object@stereo) < 2))
        return("slot 'stereo' of a Wave object must be a logical of length 1")
    if(object@stereo){
        if(length(object@left) != length(object@right))
            return("both channels of Wave objects must have the same length")
    }
    else if(length(object@right))
        return("'right' channel of a wave object is not supposed to contain data if slot stereo==FALSE")
    return(TRUE)
})

setValidity("WaveMC", 
function(object){
    if(!(mode(object@.Data) == "numeric")) return("channels of WaveMC objects must be numeric")
    return(TRUE)
})





setMethod("[", signature(x = "Wave"),
function(x, i, j, ..., drop=FALSE){
    if(!is(x, "Wave")) 
        stop("'x' needs to be of class 'Wave'")
    validObject(x)
    x@left <- x@left[i]
    if(x@stereo)
        x@right <- x@right[i]
    if(missing(j)) return(x)

    j <- gsub("1", "left", j)
    j <- gsub("2", "right", j)      
    if(!("right" %in% j))
        x <- channel(x, "left")
    if(!("left" %in% j))
        x <- channel(x, "right")
    if(length(j)==2 && j == c("right", "left"))
        x <- channel(x, "mirror")
    return(x)
})



setMethod("[", signature(x = "WaveMC"),
function(x, i, j, ..., drop=FALSE){
    if(!is(x, "WaveMC")) 
        stop("'x' needs to be of class 'WaveMC'")
    validObject(x)
    x@.Data <- x@.Data[i,j,drop=FALSE]
    return(x)
})


##########
# Wave object generating functions
setGeneric("Wave",
function(left, ...) standardGeneric("Wave"))

setGeneric("WaveMC",
function(data, ...) standardGeneric("WaveMC"))


setMethod("Wave", signature(left = "numeric"), 
function(left, right = numeric(0), samp.rate = 44100, bit = 16, pcm = TRUE, ...){
    if(missing(samp.rate)) 
        warning("'samp.rate' not specified, assuming 44100Hz")
    if(missing(bit)) 
        warning("'bit' not specified, assuming 16bit")
    return(
        new("Wave", stereo = length(right) > 0, samp.rate = samp.rate, 
            bit = bit, left = left, right = right, pcm = pcm, ...))
})

setMethod("WaveMC", signature(data = "numeric"), 
function(data = numeric(0), ...){
    WaveMC(as.matrix(data), ...)
})

setMethod("Wave", signature(left = "WaveMC"), 
function(left, ...)
    as(left, "Wave")
)

setMethod("Wave", signature(left = "matrix"), 
function(left, ...)
    Wave(as.data.frame(left), ...)
)

setMethod("WaveMC", signature(data = "matrix"), 
function(data = matrix(numeric(0), 0, 0), samp.rate = 44100, bit = 16, pcm = TRUE, ...){
    if(missing(samp.rate)) 
        warning("'samp.rate' not specified, assuming 44100Hz")
    if(missing(bit)) 
        warning("'bit' not specified, assuming 16bit")
    return(
        new("WaveMC", .Data = data, samp.rate = samp.rate, 
            bit = bit, pcm = pcm, ...))
})

setMethod("Wave", signature(left = "data.frame"), 
function(left, ...)
    Wave(as.list(left), ...)
)

setMethod("WaveMC", signature(data = "data.frame"), 
function(data, ...)
    WaveMC(as.matrix(data), ...)
)


setMethod("WaveMC", signature(data = "Wave"), 
function(data, ...)
    as(data, "WaveMC")
)


setMethod("Wave", signature(left = "list"), 
function(left, ...){
    if(length(left) > 1){
        if(length(left) > 2){    
            warning("Object has ", length(left), " channels, using only the first 2: not more than 2 channels are supported by the Wave class, use the WaveMC class instead.")
        }
        if(all(c("left", "right") %in% names(left)))
            Wave(left$left, left$right, ...)
        else 
            Wave(left[[1]], left[[2]], ...)
    }
    else Wave(left[[1]], ...)
})


setMethod("WaveMC", signature(data = "list"), 
function(data, ...){
   WaveMC(sapply(data, I), ...)
})



setAs("matrix", "Wave", function(from, to) Wave(from))
setAs("matrix", "WaveMC", function(from, to) WaveMC(from))
setAs("data.frame", "Wave", function(from, to) Wave(from))
setAs("data.frame", "WaveMC", function(from, to) WaveMC(from))
setAs("list", "Wave", function(from, to) Wave(from))
setAs("list", "WaveMC", function(from, to) WaveMC(from))
setAs("numeric", "Wave", function(from, to) Wave(from))
setAs("numeric", "WaveMC", function(from, to) WaveMC(from))

setAs("Wave", "data.frame", 
function(from, to){
    dat <- if(from@stereo) data.frame(left = from@left, right = from@right) 
           else data.frame(mono = from@left)
    return(dat)
})
setAs("WaveMC", "data.frame", 
function(from, to){
    as.data.frame(from)
})
setAs("Wave", "matrix", function(from, to) 
    return(as(as(from, "data.frame"), "matrix")))
setAs("WaveMC", "matrix", function(from, to) 
    return(from@.Data))
setAs("WaveGeneral", "list", function(from, to)
    return(as(as(from, "data.frame"), "list")))

setAs("Wave", "WaveMC", function(from, to){
    WaveMC(data = cbind(from@left, from@right), samp.rate = from@samp.rate, bit = from@bit, pcm = from@pcm)
})

setAs("WaveMC", "Wave", function(from, to){
    Wave(from@.Data, samp.rate = from@samp.rate, bit = from@bit, pcm = from@pcm)
})



setMethod("show", signature(object = "Wave"), 
function(object){
    l <- length(object@left)
    cat("\nWave Object")
    cat("\n\tNumber of Samples:     ", l)
    cat("\n\tDuration (seconds):    ",
        round(l / object@samp.rate, 2))
    cat("\n\tSamplingrate (Hertz):  ", object@samp.rate)
    cat("\n\tChannels (Mono/Stereo):",
        if(object@stereo) "Stereo" else "Mono")
    cat("\n\tPCM (integer format):  ", object@pcm)
    cat("\n\tBit (8/16/24/32/64):   ", object@bit, "\n\n")
})

setMethod("show", signature(object = "WaveMC"), 
function(object){
    l <- nrow(object)
    cat("\nWaveMC Object")
    cat("\n\tNumber of Samples:     ", l)
    cat("\n\tDuration (seconds):    ",
        round(l / object@samp.rate, 2))
    cat("\n\tSamplingrate (Hertz):  ", object@samp.rate)
    cat("\n\tNumber of channels:    ",
        ncol(object))
    cat("\n\tPCM (integer format):  ", object@pcm)
    cat("\n\tBit (8/16/24/32/64):   ", object@bit, "\n\n")
})


setMethod("summary", signature(object = "Wave"), 
function(object, ...){
    l <- length(object@left)
    cat("\nWave Object")
    cat("\n\tNumber of Samples:     ", l)
    cat("\n\tDuration (seconds):    ",
        round(l / object@samp.rate, 2))
    cat("\n\tSamplingrate (Hertz):  ", object@samp.rate)
    cat("\n\tChannels (Mono/Stereo):",
        if(object@stereo) "Stereo" else "Mono")
    cat("\n\tPCM (integer format):  ", object@pcm)
    cat("\n\tBit (8/16/24/32/64):   ", object@bit)
    cat("\n\nSummary statistics for channel(s):\n\n")
    if(object@stereo)
        print(rbind(left = summary(object@left), right = summary(object@right)))
    else print(summary(object@left))
    cat("\n\n")
})



setMethod("summary", signature(object = "WaveMC"), 
function(object, ...){
    l <- nrow(object)
    cat("\nWaveMC Object")
    cat("\n\tNumber of Samples:     ", l)
    cat("\n\tDuration (seconds):    ",
        round(l / object@samp.rate, 2))
    cat("\n\tSamplingrate (Hertz):  ", object@samp.rate)
    cat("\n\tNumber of channels:    ",
        ncol(object))
    cat("\n\tPCM (integer format):  ", object@pcm)
    cat("\n\tBit (8/16/24/32/64):   ", object@bit)
    cat("\n\nSummary statistics for channel(s):\n\n")
    print(summary(object@.Data))
    cat("\n\n")
})
