`logitreg.analog` <- function(object, groups, k = 1, ...) {
    if(is.null(object$train))
        stop("'object$train' missing. Refit 'object' with argument 'keep.train = TRUE'")
    out <- logitreg(object$train, groups = groups, k = k, ...)
    if(!is.null(attr(object, "method"))) {
        out$method <- attr(object, "method")
        ##attr(models, "method") <- attr(object, "method")
    }
    out
}
