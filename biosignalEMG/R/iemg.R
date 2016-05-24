iemg <- function(data, calliemg, resetpoints, samplingrate = 0, units = "", data.name = "") {
    if (missing(data)) 
        stop("'data' argument is not specified")
    if (!is.vector(data) | !(mode(data) == "numeric")) 
        stop("data must be a vector of numeric values")
    if (!is.numeric(samplingrate) | (samplingrate < 0)) 
        stop("Sampling rate must be 0 (for unknown sampling rate) or a positive number")
    if (!is.call(calliemg)) 
        stop("'iemg' is not intended to be called durectly, but indirectly from the 'integration' function")
    
    object <- list(values = data, call = calliemg, reset.points = resetpoints, units = units, 
        samplingrate = samplingrate, data.name = data.name)
    class(object) <- "iemg"
    return(object)
} 
