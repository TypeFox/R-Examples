onoff_singlethres <- function(data, channel, t = 0.05, data.name) {
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
    if (!is.numeric(t)) 
        stop("The threshold 't' must be a numeric value.")
    
    emgma <- envelope(data, method = "MA", wsize = 60)
    detected <- as.numeric(emgma$values > t)
    detected[is.na(detected)] <- 0
    return(detected)
} 
