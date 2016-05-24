distCalcRange <- function(densFn, param = NULL, tol = 10^(-5), ...){
    dfun <- match.fun(paste("d", densFn, sep = ""))
    mode <- distMode(densFn, param = param, ...)
    if (densFn == "skewhyp"){
        l <- list(...)
        delta <- ifelse(is.null(param), l$delta, param[2])
    }

    if(densFn != "skewhyp")
        stepHigh <- distStepSize(densFn,
                                 param = param, side = "right", ...)[1]
    xHigh <- ifelse(densFn == "skewhyp", mode + delta, mode + stepHigh)
    if (densFn == "skewhyp"){
        if(is.null(param)){
            while (dfun(xHigh, ...) > tol){
                xHigh <- xHigh +
                    distStepSize(densFn, param = param,
                                 dist = (dfun(xHigh, ...) - tol),
                                 side = "right", ...)[1]
            }
        } else {
            while (dfun(xHigh, param = param) > tol){
                xHigh <- xHigh +
                    distStepSize(densFn, param = param,
                                 dist = (dfun(xHigh, param = param) - tol),
                                 side = "right", ...)[1]
            }
        }
    } else {
        if(is.null(param)){
            while (dfun(xHigh, ...) > tol){
                xHigh <- xHigh + stepHigh
            }
        } else {
            while (dfun(xHigh, param = param) > tol){
                xHigh <- xHigh + stepHigh
            }
        }
    }

    if(densFn != "skewhyp")
        stepLow <- distStepSize(densFn,
                                param = param, side = "left", ...)[1]
    xLow <- ifelse(densFn == "skewhyp", mode - delta, mode - stepLow)
    if (densFn == "skewhyp"){
        if(is.null(param)){
            while (dfun(xLow, ...) > tol){
                xLow <- xLow -
                    distStepSize(densFn, param = param,
                                 dist = (dfun(xLow, ...) - tol),
                                 side = "left", ...)[1]
            }
        } else {
            while (dfun(xLow, param = param) > tol){
                xLow <- xLow -
                    distStepSize(densFn, param = param,
                                 dist = (dfun(xLow, param = param) - tol),
                                 side = "left", ...)[1]
            }
        }
    } else {
        if(is.null(param)){
            while (dfun(xLow, ...) > tol){
                xLow <- xLow - stepLow
            }
        } else {
            while (dfun(xLow, param = param) > tol){
                xLow <- xLow - stepLow
            }
        }
    }
    if(is.null(param)){
        zeroFun <- function(x) {
            dfun(x, ...) - tol
        }
    } else {
        zeroFun <- function(x) {
            dfun(x, param = param) - tol
        }
    }
    xUpper <- uniroot(zeroFun, c(mode,xHigh))$root
    xLower <- uniroot(zeroFun, c(xLow, mode))$root

    range <- c(xLower, xUpper)

    return(range)

}



