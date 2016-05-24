print.summary.emg <- function(x, ...) {
    object <- x
    cat("EMG Object")
    cat("\n\tTotal number of samples:", object$samples)
    cat("\n\tNumber of channels:", object$channels)
    if (object$samplingrate != 0) {
        cat("\n\tDuration (seconds): ", format(object$duration))
        cat("\n\tSamplingrate (Hertz): ", format(object$samplingrate))
    }
    if (object$channels == 1) {
        if (object$data.name != "" | object$units != "") {
            cat("\n\tChannel information:")
            if (object$data.name != "") 
                cat("\n\t", object$data.name)
            if (object$units != "") 
                cat("\n\tUnits: ", object$units)
            cat("\n")
        } else {
            cat("No additional information")
        }
    } else {
        cat("\n\tChannels information:")
        for (i in 1:object$channels) {
            cat("\n\t\t", i, ":")
            if (object$data.name[i] != "" | object$units[i] != "") {
                if (object$data.name[i] != "") 
                  cat("", object$data.name[i])
                if (object$units[i] != "") 
                  cat(", units: ", object$units[i])
                cat("\n")
            } else {
                cat("No additional information")
            }
        }
    }
} 
