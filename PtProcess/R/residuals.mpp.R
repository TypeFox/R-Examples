residuals.mpp <- function(object, ...){
    data <- object$data
    params <- object$params
    gparams <- eval(object$gmap)
    gif <- object$gif
    times <- c(0, data$time)
    tau <- rep(NA, length(times)-1)
    for (i in 2:length(times))
        tau[i-1] <- gif(data = data, evalpts = NULL, params = 
                        gparams, TT = c(times[i-1], times[i]))
    return(ts(data=cumsum(tau)))
}

