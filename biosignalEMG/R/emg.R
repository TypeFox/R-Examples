emg <- function(data, samplingrate = 0, units = "", data.name = "") {
    if (missing(data)) 
        stop("'data' argument is not specified")
    if ((!is.vector(data) & !is.matrix(data)) | !(mode(data) == "numeric")) 
        stop("data must be a vector or a matrix of numeric values")
    if (is.null(samplingrate)) 
        samplingrate <- 0
    if (!is.numeric(samplingrate) | (samplingrate < 0)) 
        stop("Sampling rate must be a positive number (or 0 for unknown sampling rate)")
    if (is.vector(data)) {
        if (length(units) != 1) 
            stop("Not valid 'units' parameter")
        if (length(data.name) != 1) 
            stop("Not valid 'data.name' parameter")
    } else {
        if (length(units) > dim(data)[2]) 
            stop("Length of 'units' parameter and number of data columns don't match")
        if (length(units) < dim(data)[2]) 
            units <- rep(units, dim(data)[2])
        if (length(data.name) > dim(data)[2]) 
            stop("Length of 'data.name' parameter and number of data columns don't match")
        if (length(data.name) < dim(data)[2]) 
            data.name <- c(data.name, rep("", dim(data)[2] - length(data.name)))
    }
    object <- list(values = data, units = units, samplingrate = samplingrate, data.name = data.name)
    class(object) <- "emg"
    return(object)
} 
