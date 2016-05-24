skewhypFit <- function (x, freq = NULL, breaks = NULL, startValues = "LA",
                        paramStart = NULL,  method = "Nelder-Mead",
                        hessian = TRUE, plots = FALSE, printOut = TRUE,
                        controlBFGS = list(maxit = 200),
                        controlNM = list(maxit = 1000),
                        maxitNLM = 1500, ...){

    xName <- paste(deparse(substitute(x), 500), collapse = "\n")

    #allow frequency data
    if (!is.null(freq)) {
        if (length(freq) != length(x)) {
            stop("vectors x and freq are not of the same length")}
        x <- rep(x, freq)
    }

    #get the starting value information
    if (startValues == "US") {
        startInfo <- skewhypFitStart(x, breaks = breaks,
                     startValues = startValues, paramStart = paramStart)
    }
    if(startValues == "LA"){
        startInfo <- skewhypFitStart(x, breaks = breaks,
                     startValues = startValues)
    }
    paramStart <- startInfo$paramStart
    svName <- startInfo$svName
    breaks <- startInfo$breaks
    empDens <- startInfo$empDens
    midpoints <- startInfo$midpoints

    #maximise likelihood
    llfunc <- function(param){
        param[2] <- exp(param[2])
        param[4] <- exp(param[4])
        -sum(dskewhyp(x, param = param, log = TRUE),
             na.rm = TRUE)}
    output <- numeric(7)
    ind <- 1:4

    if (method == "BFGS") {
        opOut <- optim(paramStart, llfunc, NULL, method = "BFGS",
            hessian = hessian, control = controlBFGS, ...)
    }
    if (method == "Nelder-Mead") {
        opOut <- optim(paramStart, llfunc, NULL, method = "Nelder-Mead",
            hessian = hessian, control = controlNM, ...)
    }
    if (method == "nlm") {
        ind <- c(2, 1, 5, 4)
        opOut <- nlm(llfunc, paramStart, hessian = hessian, iterlim = maxitNLM,
            ...)
    }

    #results
    param <- as.numeric(opOut[[ind[1]]])[1:4]
    param[2] <- exp(param[2])
    param[4] <- exp(param[4])
    names(param) <- c("mu", "delta", "beta", "nu")
    maxLik <- -(as.numeric(opOut[[ind[2]]]))
    conv <- as.numeric(opOut[[ind[4]]])
    iter <- as.numeric(opOut[[ind[3]]])[1]
    paramStart <- c(paramStart[1], exp(paramStart[2]), paramStart[3],
        exp(paramStart[4]))#unlog the logged parameters
    fitResults <- list(param=param , maxLik = maxLik,
        hessian = if (hessian) opOut$hessian else NULL,
        method = method, conv = conv, iter = iter, x=x, xName = xName,
        paramStart = paramStart, svName = svName, startValues = startValues,
        breaks = breaks, midpoints = midpoints, empDens = empDens)

    class(fitResults) <- "skewhypFit"
    if (printOut == TRUE) {
        print.skewhypFit(fitResults, ...)
    }
    if (plots == TRUE) {
        plot.skewhypFit(fitResults, ...)
    }
    return(fitResults)
}


