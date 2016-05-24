setGeneric("noSilence",
function(object, zero = 0, level = 0, where = c("both", "start", "end")) 
    standardGeneric("noSilence"))


setMethod("noSilence", signature(object = "Wave"),
function(object, zero = 0, level = 0, where = c("both", "start", "end")){
    validObject(object)
    where <- match.arg(where)
    starts <- 1
    ends <- length(object)    
    
    Lzero <- if(is.na(zero)) mean(object@left) else zero
    Rzero <- if(object@stereo && is.na(zero)) mean(object@right) else zero
    if(!is.numeric(Lzero) || length(Lzero) != 1)
        stop("zero must be either NA or a numeric value of length 1")
    if(!is.numeric(level) || length(level) != 1)
        stop("level must be a numeric value of length 1")
    NoSilence <- which(abs(object@left - Lzero) > level)
    if(length(NoSilence))
        NoSilence <- range(NoSilence) 
    if(object@stereo){
        temp <- which(abs(object@right - Rzero) > level)
        if(length(temp))
            NoSilence <- range(NoSilence, temp) 
    }

    if(length(NoSilence) > 1){
        if(where %in% c("both", "start"))
            starts <- NoSilence[1]
        if(where %in% c("both", "end"))
            ends <- NoSilence[2]
        return(object[starts:ends])
    } else{ 
        warning("returned object is empty (was completely silent)")
        return(object[NoSilence])
    }
})





setMethod("noSilence", signature(object = "WaveMC"),
function(object, zero = 0, level = 0, where = c("both", "start", "end")){
    validObject(object)
    where <- match.arg(where)
    starts <- 1
    ends <- length(object)    
    
    Czero <- if(is.na(zero)) colMeans(object) else rep(zero, ncol(object))
    if(!is.numeric(Czero))
        stop("zero must be either NA or numeric")
    if(!is.numeric(level))
        stop("level must be a numeric")
    NoSilence <- range(which(apply(abs(t(object) - Czero) > level, 2, any)))

    if(length(NoSilence) > 1){
        if(where %in% c("both", "start"))
            starts <- NoSilence[1]
        if(where %in% c("both", "end"))
            ends <- NoSilence[2]
        return(object[starts:ends])
    } else{ 
        warning("returned object is empty (was completely silent)")
        return(object[NoSilence])
    }
})
