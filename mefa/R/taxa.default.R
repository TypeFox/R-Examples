taxa.default <- function(x, summaries=FALSE, ...) {
    if (inherits(x, "Mefa")) return(x@xtab)
    if (inherits(x, "mefa")) {
        if (is.null(x$taxa))
            return(NULL) else if (summaries)
                return(as.data.frame(x, fun=mts, ...)) else return(x$taxa)
    }
    stop("not mefa class")
}
