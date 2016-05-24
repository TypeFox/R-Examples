pDist <- function(densFn = "norm", q, param = NULL,
                  subdivisions = 100, lower.tail = TRUE,
                  intTol = .Machine$double.eps^0.25,
                  valueOnly = TRUE, ...)
{
    CALL <- match.call()
    dfun <- match.fun(paste("d", densFn, sep = ""))
    mode <- distMode(densFn, param = param, ...)
    ## match the density function for different distributions
    qLess <- which((q <= mode)&(is.finite(q)))
    ## when q is less than mode
    qGreater <- which((q > mode)&(is.finite(q)))
    ## when q is greater than mode
    prob <- rep(NA, length(q))
    err <- rep(NA, length(q))
    prob[q == -Inf] <- 0
    prob[q == Inf] <- 0
    err[q %in% c(-Inf, Inf)] <- 0

    ## Integrate density from Inf / -Inf to the mode
    for (i in qLess){
        if(is.null(param)){
            intRes <- integrate(dfun, -Inf, q[i],
                                subdivisions = subdivisions,
                                rel.tol = intTol, ...)
        } else {
            intRes <- integrate(dfun, -Inf, q[i], param = param,
                                subdivisions = subdivisions,
                                rel.tol = intTol, ...)
        }
        prob[i] <- intRes$value
        err[i] <- intRes$abs.error
    }

    for (i in qGreater){
        if(is.null(param)){
            intRes <- integrate(dfun, q[i], Inf,
                                subdivisions = subdivisions,
                                rel.tol = intTol, ...)
        } else {
            intRes <- integrate(dfun, q[i], Inf, param = param,
                                subdivisions = subdivisions,
                                rel.tol = intTol, ...)
        }
        prob[i] <- intRes$value
        err[i] <- intRes$abs.error
    }

    if (lower.tail == TRUE)
   {
       prob[q > mode] <- 1 - prob[q > mode]
   }
   else
   {
       prob[q <= mode] <- 1 - prob[q <= mode]
   }

    # Return Value:
    ifelse(valueOnly, return(prob),
           return(list(value = prob, error = err)))

}


qDist <- function(densFn = "norm", p, param = NULL,
                  lower.tail = TRUE, method = "spline",  nInterpol = 501,
                  uniTol = .Machine$double.eps^0.25, subdivisions = 100,
                  intTol = uniTol, ...)
{
    CALL <- match.call()
    if (!lower.tail) {
        p <- 1 - p
        lower.tail <- TRUE
    }
    mode <- distMode(densFn, param = param, ...)
    pMode <- pDist(densFn, q = mode, param = param, intTol = intTol, ...)
    quant <- rep(NA, length(p))
    invalid <- which((p < 0) | (p > 1))
    pFinite <- which((p > 0) & (p < 1))
    if(densFn == "skewhyp"){
        l <- list(...)
        delta <- ifelse(is.null(param), l$delta, param[2])
    } else delta <- 0

    xRange <- distCalcRange(densFn, param = param, tol = 10^(-5), ...)

    if (method == "integrate"){
        less <- which((p <= pMode) & (p >.Machine$double.eps^7.5))
        quant <- ifelse(p <= .Machine$double.eps^5, -Inf, quant)
        if (length(less) > 0){
            pLow <- min(p[less])
            step <- distStepSize(densFn, param = param,
                                 dist = delta, side = "left", ...)[1]
            xLow <- mode - step
            if (densFn == "skewhyp"){
                while(pDist(densFn, xLow,
                            param = param, intTol = intTol, ...) >= pLow)
                {
                    xLow <- xLow -
                        distStepSize(densFn, param = param,
                                     dist = delta, side = "left", ...)[1]
                }
            } else {
                while(pDist(densFn, xLow,
                            param = param, intTol = intTol, ...) >= pLow)
                {
                    xLow <- xLow - step
                }
            }
            xRange <- c(xLow, mode)
            zeroFn <- function(x, p){
                return(pDist(densFn, x, param = param,
                             intTol = intTol, ...) - p)
            }
            for (i in less)
            {
                quant[i] <- uniroot(zeroFn, p = p[i],
                                    interval = xRange, tol = uniTol)$root
            }
        }

        greater <- which ((p > pMode) & (p < (1 - .Machine$double.eps^7.5)))
        p[greater] <- 1 - p[greater]
        quant <- ifelse(p >= (1 - .Machine$double.eps^5), Inf, quant)
        if (length(greater) > 0){
            pHigh <- min(p[greater])
            step <-  distStepSize(densFn, param = param,
                                  dist = delta, side = "right", ...)[1]
            xHigh <- mode + step
            if (densFn == "skewhyp"){
                while (pDist(densFn, xHigh, param = param, intTol = intTol,
                             lower.tail = FALSE, ...) >= pHigh)
                {
                    xHigh <- xHigh +
                        distStepSize(densFn, param = param,
                                     dist = delta, side = "left", ...)[1]
                }
            } else {
                while (pDist(densFn, xHigh, param = param, intTol = intTol,
                             lower.tail = FALSE, ...) >= pHigh)
                {
                    xHigh <- xHigh + step
                }
            }
            xRange <- c(mode, xHigh)
            zeroFn <- function(x, p){
                return(pDist(densFn, x, param = param, intTol = intTol,
                             lower.tail = FALSE, ...) - p)
            }
            for (i in greater)
            {
                quant[i] <- uniroot(zeroFn, p = p[i],
                                    interval = xRange, tol = uniTol)$root
            }
        }
    } else if (method == "spline"){
        inRange <- which((p > pDist(densFn, xRange[1],
                                    param = param, intTol = intTol, ...)) &
                         (p < pDist(densFn, xRange[2],
                                    param = param, intTol = intTol, ...)))
        small <- which((p <= pDist(densFn, xRange[1],
                        param = param, intTol = intTol, ...)) & (p >= 0))
        large <- which((p >= pDist(densFn, xRange[2],
                        param = param, intTol = intTol, ...)) & (p <= 1))
        extreme <- c(small, large)
        xVals <- seq(xRange[1], xRange[2], length.out = nInterpol)
        yVals <- pDist(densFn, xVals, param = param,
                       subdivisions = subdivisions, intTol = intTol, ...)
        splineFit <- splinefun(xVals, yVals)
        zeroFn <- function(x, p){
            return(splineFit(x) - p)
        }

        for (i in inRange){
            quant[i] <- uniroot(zeroFn, p = p[i],
                                interval = xRange, tol = uniTol)$root
        }

        if (length(extreme) > 0){
            quant[extreme] <- qDist(densFn, p[extreme], param = param,
                                    lower.tail = lower.tail,
                                    method = "integrate",
                                    nInterpol = nInterpol, uniTol = uniTol,
                                    subdivisions = subdivisions,
                                    intTol = intTol, ...)
        }
    }


    return(quant)
}



