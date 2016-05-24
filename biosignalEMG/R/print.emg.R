print.emg <- function(x, ...) {
    object <- x
    cat("EMG Object")
    if (object$data.name != "") 
        cat(" : ", object$data.name)
    cat("\n\tNumber of samples:", length(object$values))
    if (object$units != "") 
        cat("\n\tUnits: ", object$units)
    if (object$samplingrate != 0) {
        cat("\n\tDuration (seconds): ", format(length(object$values)/object$samplingrate))
        cat("\n\tSamplingrate (Hertz): ", format(object$samplingrate))
    }
    cat("\nValues:\n")
    print(object$values)
} 
