envelope <- function(data, channel, method = c("MA", "RMS"), wsize, data.name, ...) {
    if (missing(data)) 
        stop("'data' argument is not specified")
    if (!is.emg(data)) 
        stop("an object of class 'emg' is required")
    method <- match.arg(method)
    if (missing(wsize)) 
        stop("Window size argument is requiered")
    args <- list(...)
    namesargs <- names(args)
    if ((method == "MA") & (!("rtype" %in% namesargs))) 
        rtype <- "fullwave" else rtype <- args$rtype
    if (!("units" %in% namesargs)) 
        units <- "samples" else units <- args$units
    if (missing(channel)) {
        if (missing(data.name)) 
            data <- extractchannel(data) else data <- extractchannel(data, data.name = data.name)
    } else {
        if (missing(data.name)) 
            data <- extractchannel(data, channel) else data <- extractchannel(data, channel, data.name)
    }
    
    if (method == "MA") {
        rsvalues <- rectification(data, rtype = rtype)
    } else {
        rsvalues <- emg((data$values - mean(data$values))^2)
    }
    evalues <- movingaverage(rsvalues, wsize = wsize, units = units)$values
    object <- emg(evalues, data$samplingrate, data$units, data$data.name)
    return(object)
} 
