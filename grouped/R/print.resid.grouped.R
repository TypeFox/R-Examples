"print.resid.grouped" <-
function(x, ...){
    if(!inherits(x, "resid.grouped"))
        stop("Use only with 'resid.grouped' objects.\n")
    res <- x$residuals
    names(res) <- x$nam.res
    print(res)
    invisible(x)
}

