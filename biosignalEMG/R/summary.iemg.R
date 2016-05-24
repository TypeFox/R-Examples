summary.iemg <- function(object, ...) {
    res <- list(data.name = object$data.name, call = object$call, samples = length(object$values), 
        units = object$units, reset.points = object$reset.points)
    if (object$samplingrate != 0) {
        res$duration <- length(object$values)/object$samplingrate
        res$samplingrate = object$samplingrate
    }
    class(res) <- "summary.iemg"
    return(res)
} 
