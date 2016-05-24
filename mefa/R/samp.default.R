samp.default <- function(x, summaries=FALSE, ...) {
    if (inherits(x, "Mefa")) return(x@xtab)
    if (inherits(x, "mefa")) {
        if (is.null(x$samp))
            return(NULL) else if (summaries)
                return(as.data.frame(x, fun=mss, ...)) else return(x$samp)
    }
    stop("not mefa class")
}

