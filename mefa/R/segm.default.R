segm.default <- function(x, segments=1:dim(x)[3], ...) {
    if (inherits(x, "Mefa")) return(x@xtab)
    if (inherits(x, "mefa")) {
        if (!all(segments %in% 1:dim(x)[3]))
            stop("'segments' out of range")
        if (is.null(x$segm))
            return(x$xtab) else return(x$segm[segments])
    }
    stop("not mefa class")
}
