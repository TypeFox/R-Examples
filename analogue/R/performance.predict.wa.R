`performance.predict.wa` <- function(object, ...) {
    if(is.null(object$performance))
        stop("No permformance statistics in 'object'")
    retval <- with(object$performance,
                   c(rmsep, r.squared, avg.bias, max.bias))
    names(retval) <- c("RMSEP","R2","Avg.Bias","Max.Bias")
    class(retval) <- "performance"
    attr(retval, "CV.method") <- object$CV.method
    retval
}
