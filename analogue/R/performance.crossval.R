`performance.crossval` <- function(object, ...) {
    retval <- object$performance
    class(retval) <- c("performance","data.frame")
    retval
}
