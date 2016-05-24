as.emg.matrix <- function(x, ...) {
    args <- list(...)
    namesargs <- names(args)
    if ("samplingrate" %in% namesargs) 
        samplingrate <- args$samplingrate else samplingrate <- 0
    if ("units" %in% namesargs) 
        units <- args$units else units <- ""
    data <- x
    if (dim(data)[2] > 1) {
        dimnames(data) <- NULL
    } else {
        attributes(data) <- NULL
    }
    object <- emg(data, samplingrate = samplingrate, units = units, data.name = colnames(x))
    return(object)
} 
