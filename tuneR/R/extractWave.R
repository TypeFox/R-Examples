extractWave <-
function(object, from = 1, to = length(object), 
    interact = interactive(), xunit = c("samples", "time"), ...){

    if(!(is(object, "Wave") || is(object, "WaveMC"))) 
        stop("'object' needs to be of class 'Wave' or 'WaveMC'")
    validObject(object)

    xunit <- match.arg(xunit)

    if(interact){
        mf <- missing(from)
        mt <- missing(to)
        if(mf || mt)
            plot(object, xunit = xunit, ...)
        if(mf){
            cat("Click for 'from'\t")
            if(.Platform$OS.type == "windows") flush.console()
            from <- locator(1)$x
            cat(from, "\n")
            abline(v = from, ...)
        }
        if(mt){
            cat("Click for 'to'  \t")
            if(.Platform$OS.type == "windows") flush.console()
            to <- locator(1)$x
            cat(to, "\n")
            abline(v = to, ...)
        }
    }
    
    if(xunit == "time"){
        to <- to * object@samp.rate
        from <- from * object@samp.rate
    }

    lo <- length(object)
    from <- max(from, 1)
    to <- min(to, lo)

    if(from > to){
        warning("'from' > 'to', object is unchanged")
        return(object)
    }
    
    return(object[seq(from, to)])
}
