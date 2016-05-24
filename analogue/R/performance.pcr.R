`performance.pcr` <- function(object, ...) {
    retval <- object$performance
    class(retval) <- c("performance","data.frame")
    retval
}
