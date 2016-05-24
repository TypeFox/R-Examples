extractchannel <- function(data, channel, data.name) {
    if (missing(data)) 
        stop("'data' argument is not specified")
    if (!is.emg(data)) 
        stop("an object of class 'emg' is required")
    if (is.vector(data$values)) {
        if (!missing(channel)) 
            if (channel != 1) 
                warning("Unused parameter 'channel'")
        values <- data$values
        if (missing(data.name)) 
            data.name <- data$data.name
        units <- data$units
    } else {
        if (missing(channel)) 
            stop("Parameter 'channel' is missing")
        v <- .validchannel(data, channel)
        if (v$val == 1) {
            values <- data$values[, v$channel]
            if (missing(data.name)) 
                data.name <- data$data.name[v$channel]
            units <- data$units[v$channel]
        } else {
            stop(v$message)
        }
    }
    object <- emg(values, data$samplingrate, units, data.name = data.name)
    return(object)
} 
