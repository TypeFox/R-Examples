setMethod("plot", signature(x = "Wspec", y = "missing"),
function(x, which = 1, type = "h", xlab = "Frequency", ylab = NULL, log = "", ...){

    if(is.null(ylab)){
        ylab <- if(x@normalize) "normalized periodogram" else "periodogram"
        if(!missing(log) && (log == "y")) ylab <- paste("log(", ylab, ")", sep = "")
    }        
    spec <- x@spec[[which]]
    plot(x@freq, spec, type = type, 
        xlab = xlab, ylab = ylab, log = log, ...)
})


setMethod("plot", signature(x = "WspecMat", y = "missing"),
function(x, xlab = "time", ylab = "Frequency", xunit = c("samples", "time"), log = "", ...){
    if(log == "z"){ 
        x@spec <- log(x@spec)
        log <- ""
    }
    xunit <- match.arg(xunit)
    if(xunit == "time"){
        x@starts <- x@starts / x@samp.rate
    }
    image(x@starts, x@freq, x@spec, xlab = xlab, ylab = ylab, log = log, ...)
})

setMethod("image", signature(x = "Wspec"),
function(x, xlab = "time", ylab = "Frequency", xunit = c("samples", "time"), log = "", ...){
    x <- as(x, "WspecMat")
    xunit <- match.arg(xunit)    
    plot(x, xlab = xlab, ylab = ylab, xunit = xunit, log = log, ...)
})
