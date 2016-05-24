print.summary.iemg <- function(x, ...) {
    object <- x
    cat("Call:\n")
    print(object$call)
    cat("\n(Integrated) EMG Object")
    if (object$data.name != "") 
        cat(" : ", object$data.name)
    cat("\n\tNumber of samples:", object$samples)
    if (!is.null(object$reset.points)) {
        cat("\n\tIntegrated signal reset at sample points:\n\t\t")
        cat(object$reset.points)
    }
    if (object$units != "") 
        cat("\n\tUnits: ", object$units)
    if (!is.null(object$samplingrate)) {
        cat("\n\tDuration (seconds): ", format(object$duration))
        cat("\n\tSamplingrate (Hertz): ", format(object$samplingrate))
    }
    cat("\n")
} 
