summary.emg <- function(object, ...) {
    res <- list(data.name = object$data.name, samples = length(object$values), units = object$units)
    
    if (is.vector(object$values)) {
        if (object$samplingrate != 0) 
            res$duration <- length(object$values)/object$samplingrate
        res$channels <- 1
    } else {
        if (object$samplingrate != 0) 
            res$duration <- dim(object$values)[1]/object$samplingrate
        res$channels <- dim(object$values)[2]
    }
    
    res$samplingrate = object$samplingrate
    class(res) <- "summary.emg"
    return(res)
} 
