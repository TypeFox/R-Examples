print.iemg <- function(x, ...) {
    object <- x
    cat("(Integrated) EMG Object")
    if (object$data.name != "") 
        cat(" : ", object$data.name)
    cat("\n\tNumber of samples:", length(object$values))
    if (!is.null(object$reset.points)) {
        cat("\n\tIntegrated signal reset at sample points:\n\t\t")
        cat(object$reset.points)
    }
    if (object$units != "") 
        cat("\n\tUnits: ", object$units)
    if (object$samplingrate != 0) {
        cat("\n\tDuration (seconds): ", format(length(object$values)/object$samplingrate))
        cat("\n\tSamplingrate (Hertz): ", format(object$samplingrate))
    }
    cat("\nValues:\n")
    print(object$values)
} 
