rectification <- function(data, channel, rtype = c("fullwave", "halfwave"), data.name, 
    ...) {
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
    rtype <- match.arg(rtype)
    
    tdata <- data$values
    rectdata <- abs(tdata)
    if (rtype == "halfwave") 
        rectdata[tdata < 0] <- 0
    object <- emg(rectdata, data$samplingrate, data$units, data$data.name)
    return(object)
} 
