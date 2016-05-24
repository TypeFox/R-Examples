xtab.default <- function(x, ...) {
    if (inherits(x, "mefa")) return(x$xtab)
    if (inherits(x, "Mefa")) return(x@xtab)
    stop("not mefa class")
}
