highpass <- function(data, channel, cutoff = 50, data.name) {
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
    
    if (data$samplingrate == 0) 
        stop("The sampling rate is requiered")
    bf <- butter(5, 2 * cutoff/data$samplingrate, type = "high")
    b <- signal::filter(bf, data$values)
    attributes(b) <- NULL
    object <- emg(b, data$samplingrate, data$units, data$data.name)
    return(object)
} 
