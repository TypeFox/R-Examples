`print.performance` <- function(x, digits = min(getOption("digits"), 4),
                                ...) {
    CV.method <- attr(x, "CV.method")
    if(inherits(x, "data.frame")) {
        print.data.frame(x, digits = digits, ...)
    } else {
        x <- zapsmall(x, digits = digits)
        perf.names <- names(x)
        attributes(x) <- NULL
        names(x) <- perf.names
        print.default(x, digits = digits, ...)
    }
    invisible(x)
}
