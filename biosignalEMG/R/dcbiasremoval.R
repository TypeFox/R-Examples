dcbiasremoval <- function(data, channel, baseline, data.name) {
    if (missing(data)) 
        stop("'data' argument is not specified")
    if (!is.emg(data)) 
        stop("an object of class 'emg' is required")
    if (missing(channel)) {
        if (missing(data.name)) 
            data <- extractchannel(data) else data <- extractchannel(data, data.name = data.name)
    } else {
        if (missing(data.name)) 
            data <- extractchannel(data, channel) else data <- extractchannel(data, channel, data.name)
    }
    if (missing(baseline)) 
        baseline <- mean(data$values)
    values <- data$values - baseline
    object <- emg(values, data$samplingrate, data$units, data$data.name)
    return(object)
} 
