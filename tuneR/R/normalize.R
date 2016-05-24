setGeneric("normalize",
  function(object, unit = c("1", "8", "16", "24", "32", "64", "0"), center = TRUE, level = 1, rescale = TRUE, pcm = object@pcm)
    standardGeneric("normalize"))

setMethod("normalize", signature(object = "Wave"),
  function(object, unit = c("1", "8", "16", "24", "32", "64", "0"), center = TRUE, level = 1, rescale = TRUE, pcm = object@pcm){   
    validObject(object)
    unit <- match.arg(unit)
    if(!(unit %in% c("1", "8", "16", "24", "32", "64", "0")))
        stop("'unit' must be either 1 (real valued norm.), 8 (norm. to 8-bit), 16 (norm. to 16-bit), 24 (...), 32 (integer or real valued norm., depends on pcm), or  64 (real valued norm.)")
    if(unit == 64 && pcm){
        warning("pcm set to FALSE since unit=64")
        pcm <- FALSE
    }
    if(unit %in% c(8, 16, 24) && !pcm){        
        warning("pcm set to TRUE since unit was one of 8, 16, or 24")
        pcm <- TRUE
    }
    if(center){
        object@left <- object@left - mean(object@left)
        object@right <- object@right - mean(object@right)    
    }
    if(object@bit == 8 && all(object@right >= 0) && all(object@left >= 0)){
        object@left <- object@left - 127
        object@right <- object@right - 127    
    }   
    if(unit != "0"){
        if(rescale){
            m <- max(abs(c(range(object@left), if(object@stereo) range(object@right))))
        } else {
            if(object@pcm) {
                m <- switch(as.character(object@bit),
                    "1" = 1,
                    "8" = 128,
                    "16" = 32768,
                    "24" = 8388608,
                    "32" = 2147483648) 
            } else m <- 1
        }
        if(!isTRUE(all.equal(m, 0))){
            object@left <- level * object@left / m
            object@right <- level * object@right / m
        }
        if(pcm){
            if(unit == "8"){
                object@left <- round(object@left * 127 + 127)
                object@right <- round(object@right * 127 + 127)
            }
            else if(unit == "16"){
                object@left <- round(object@left * 32767)
                object@right <- round(object@right * 32767)
            }    
            else if(unit == "24"){
                object@left <- round(object@left * 8388607)
                object@right <- round(object@right * 8388607)
            }    
            else if(unit == "32"){
                object@left <- round(object@left * 2147483647)
                object@right <- round(object@right * 2147483647)
            }    
        }
    }
    unit <- as.integer(unit)
    if(unit){
        object@bit <- if(unit==1) 32 else unit
        object@pcm <- pcm
    }
    return(object)
})





setMethod("normalize", signature(object = "WaveMC"),
  function(object, unit = c("1", "8", "16", "24", "32", "64", "0"), center = TRUE, level = 1, rescale = TRUE, pcm = object@pcm){   
    validObject(object)
    unit <- match.arg(unit)
    if(!(unit %in% c("1", "8", "16", "24", "32", "64", "0")))
        stop("'unit' must be either 1 (real valued norm.), 8 (norm. to 8-bit), 16 (norm. to 16-bit), 24 (...), 32 (integer or real valued norm., depends on pcm), or  64 (real valued norm.)")
    if(unit == 64 && pcm){
        warning("pcm set to FALSE since unit=64")
        pcm <- FALSE
    }
    if(unit %in% c(8, 16, 24) && !pcm){        
        warning("pcm set to TRUE since unit was one of 8, 16, or 24")
        pcm <- TRUE
    }
    if(center){
        object@.Data <- scale(object@.Data, scale=FALSE, center=TRUE)
    }
    if(object@bit == 8 && all(object >= 0)){
        object <- object - 127
    }
    if(unit != "0"){
        if(rescale){
            m <- max(abs(range(object)))
        } else {
            if(object@pcm) {
                m <- switch(as.character(object@bit),
                    "1" = 1,
                    "8" = 128,
                    "16" = 32768,
                    "24" = 8388608,
                    "32" = 2147483648) 
            } else m <- 1
        }
        if(!isTRUE(all.equal(m, 0))){
            object <- level * object / m
        }
        if(pcm && unit %in% c("8", "16", "24", "32")){
            object <- round(
                switch(as.character(unit),
                    "8" = {object * 127 + 127},
                    "16" = {object * 32767},
                    "24" = {object * 8388607},
                    "32" = {object * 2147483647}))
        }
    }
    unit <- as.integer(unit)
    if(unit){
        object@bit <- if(unit==1) 32 else unit
        object@pcm <- pcm
    }
    return(object)
})
