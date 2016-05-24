as.emg.numeric <- function(x, ...) {
    args <- list(...)
    namesargs <- names(args)
    if ("samplingrate" %in% namesargs) 
        samplingrate <- args$samplingrate else samplingrate <- 0
    if ("units" %in% namesargs) 
        units <- args$units else units <- ""
    data <- x
    attributes(data) <- NULL
    object <- emg(data, samplingrate = samplingrate, units = units, data.name = deparse(substitute(x)))
    return(object)
} 
