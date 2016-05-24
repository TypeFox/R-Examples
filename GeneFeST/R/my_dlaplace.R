my_dlaplace <-function (x, location = 0, scale = 1, log = FALSE) 
{
    if (!is.logical(log.arg <- log) || length(log) != 1) 
        stop("bad input for argument 'log'")
    rm(log)
    logdensity <- (-abs(x - location)/scale) - log(2 * scale)
    if (log.arg) 
        logdensity
    else exp(logdensity)
}

