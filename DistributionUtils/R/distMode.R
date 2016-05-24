distMode <- function(densFn, param = NULL, ...){
    dfun <- match.fun(paste("d", densFn, sep = ""))
    if(densFn == "skewhyp"){
        l <- list(...)
        delta <- ifelse(is.null(param), l$delta, param[2])
    } else delta <- 0
    median <- distStepSize(densFn, param = param,
                           dist = delta, side = "left", ...)[2]
    if(is.null(param)){
        modefun <- function(x){
            log(dfun(x, ...))
        }
    } else {
        modefun <- function(x){
            log(dfun(x, param = param))
        }
    }
    stepHigh <- distStepSize(densFn, param = param,
                             dist = delta, side = "right", ...)[1]
    xHigh <- median + stepHigh
    if (densFn == "skewhyp"){
        if(is.null(param)){
            while(dfun(xHigh, ...) > dfun(median, ...)){
                xHigh <- xHigh +
                    distStepSize(densFn, param = param,
                                 dist = delta, side = "right", ...)[1]
            }
        } else {
            while(dfun(xHigh, param = param) > dfun(median, param = param)){
                xHigh <- xHigh +
                    distStepSize(densFn, param = param,
                                 dist = delta, side = "right", ...)[1]
            }
        }
    } else {
        if(is.null(param)){
            while(dfun(xHigh, ...) > dfun(median, ...)){
                xHigh <- xHigh + stepHigh
            }
        } else {
            while(dfun(xHigh, param = param) > dfun(median, param = param)){
                xHigh <- xHigh + stepHigh
            }
        }
    }
    stepLow <-  distStepSize(densFn, param = param,
                             dist = delta, side = "left", ...)[1]
    xLow <- median - stepLow
    if (densFn == "skewhyp"){
        if(is.null(param)){
            while(dfun(xLow, ...) > dfun(median, ...)){
                xLow <- xLow -
                    distStepSize(densFn, param = param,
                                 dist = delta, side = "left", ...)[1]
            }
        } else {
            while(dfun(xLow, param = param) > dfun(median, param = param)){
                xLow <- xLow -
                    distStepSize(densFn, param = param,
                                 dist = delta, side = "left", ...)[1]
            }
        }
    } else {
        if(is.null(param)){
            while(dfun(xLow, ...) > dfun(median, ...)){
                xLow <- xLow - stepLow
            }
        } else {
            while(dfun(xLow, param = param) > dfun(median, param = param)){
                xLow <- xLow - stepLow
            }
        }
    }
    range <- c(xLow, xHigh)
    optResult <- optimize(f = modefun, interval = range, maximum = TRUE)
    mode <- optResult$maximum
    mode
}




