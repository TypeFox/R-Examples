residuals.linksrm <- function(object, ...){
    data <- object$data
    params <- object$params
    gparams <- eval(object$gmap)
    gif <- object$gif
    times <- c(0, data$time)
    tau <- NULL
    for (i in 2:length(times))
        tau <- rbind(tau, gif(data = data, evalpts = NULL,
                params = gparams, TT = c(times[i-1], times[i])))
    #   tau is a matrix with ncol being number of regions
    #   and nrow being total number of events
    x <- list()
    for (region in 1:max(data[,"region"])){
        tau[,region] <- cumsum(tau[,region])
        x[[region]] <- ts(tau[(data[,"region"]==region), region])
        names(x[[region]]) <- NULL
    }
    return(x)
}

