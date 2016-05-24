require(methods)
##########
# define class Wspec
setClass("Wspec",
    representation = representation(
        freq = "numeric",
        spec = "list", 
        # coh, phase,
        kernel = "ANY",
        df = "numeric",
        taper = "numeric",
        width = "numeric",
        overlap = "numeric",
        normalize = "logical",
        starts = "numeric",
        stereo = "logical",
        samp.rate = "numeric",
        variance = "numeric",
        energy = "numeric"),
    prototype = prototype(taper = 0, 
        stereo = FALSE, samp.rate = 44100))


setMethod("[", signature(x = "Wspec"),
function(x, i, j, ..., drop=FALSE){
    if(!is(x, "Wspec")) 
        stop("'x' needs to be of class 'Wspec'")
    validObject(x)

    x@spec <- x@spec[i]
    x@starts <- x@starts[i]
    x@variance <- x@variance[i]
    x@energy <- x@energy[i]
    return(x)
})


setMethod("show", signature(object = "Wspec"), 
function(object){
    l <- length(object@freq)
    cat("Wspec Object (use summary() for more details)\n\n")
    cat("Number of Periodograms:", length(object@spec), "\n")
    cat("Estimated at", l, "Frequencies:", 
        object@freq[1], "...", object@freq[l], "\n\n")
    cat("Further parameters:\n")
    cat("width:  ", object@width,   "\n")
    cat("overlap:", object@overlap, "\n")
    cat("normal.:", object@normalize, "\n\n")
})


setMethod("summary", signature(object = "Wspec"),
function(object, ...){
    l <- length(object@freq)
    cat("Wspec Object", 
        if(!is.null(object@kernel)) 
            "(see object@kernel for details on the kernel)", 
        "\n\n")
    cat("Number of Periodograms:", length(object@spec), "\n")
    cat("Estimated at", l, "Frequencies:", 
        object@freq[1], "...", object@freq[l], "\n\n")
    cat("Further parameters:\n")
    cat("df:     ", object@df,      "\n")
    cat("taper:  ", object@taper,   "\n")
    cat("width:  ", object@width,   "\n")
    cat("overlap:", object@overlap, "\n")
    cat("normal.:", object@normalize, "\n\n")

    cat("Properties of the Wave object:\n")
    cat("Sampling rate:", object@samp.rate, "\n\n")
})



##########
# define class Wspec
setClass("WspecMat",
    representation = representation(
        freq = "numeric",
        spec = "matrix",
        # coh, phase,
        kernel = "ANY",
        df = "numeric",
        taper = "numeric",
        width = "numeric",
        overlap = "numeric",
        normalize = "logical",
        starts = "numeric",
        stereo = "logical",
        samp.rate = "numeric",
        variance = "numeric",
        energy = "numeric"),
    prototype = prototype(taper = 0, 
        stereo = FALSE, samp.rate = 44100))



setMethod("[", signature(x = "WspecMat"),
function(x, i, j, ..., drop=FALSE){
    if(!is(x, "WspecMat")) 
        stop("'x' needs to be of class 'WspecMat'")
    validObject(x)

    x@spec <- x@spec[i,j]
    x@starts <- x@starts[i]
    x@variance <- x@variance[i]
    x@energy <- x@energy[i]
    return(x)
})


setMethod("show", signature(object = "WspecMat"), 
function(object){
    l <- length(object@freq)
    cat("WspecMat Object (use summary() for more details)\n\n")
    cat("Number of Periodograms:", nrow(object@spec), "\n")
    cat("Estimated at", l, "Frequencies:", 
        object@freq[1], "...", object@freq[l], "\n\n")
    cat("Further parameters:\n")
    cat("width:  ", object@width,   "\n")
    cat("overlap:", object@overlap, "\n")
    cat("normal.:", object@normalize, "\n\n")
})


setMethod("summary", signature(object = "WspecMat"),
function(object, ...){
    l <- length(object@freq)
    cat("WspecMat Object", 
        if(!is.null(object@kernel)) 
            "(see object@kernel for details on the kernel)", 
        "\n\n")
    cat("Number of Periodograms:", nrow(object@spec), "\n")
    cat("Estimated at", l, "Frequencies:", 
        object@freq[1], "...", object@freq[l], "\n\n")
    cat("Further parameters:\n")
    cat("df:     ", object@df,      "\n")
    cat("taper:  ", object@taper,   "\n")
    cat("width:  ", object@width,   "\n")
    cat("overlap:", object@overlap, "\n")
    cat("normal.:", object@normalize, "\n\n")

    cat("Properties of the Wave object:\n")
    cat("Sampling rate:", object@samp.rate, "\n\n")
})

setAs("Wspec", "WspecMat", function(from, to) {
    spec <- from@spec
    l <- sapply(spec, length)
    ml <- max(l)
    for(i in which(l < ml))
        spec[[i]] <- rep(spec[[i]], length.out = ml)
    spec <- matrix(unlist(spec), nrow = length(spec), byrow=TRUE)
    return(new("WspecMat", 
        freq = from@freq,
        spec = spec, 
        # coh, phase,
        kernel = from@kernel,
        df = from@df,
        taper = from@taper,
        width = from@width,
        overlap = from@overlap,
        normalize = from@normalize,
        starts = from@starts,
        stereo = from@stereo,
        samp.rate = from@samp.rate,
        variance = from@variance,
        energy = from@energy))
})
