####### QQ Plot function ##################################################
qqskewhyp <- function(y, mu = 0, delta = 1, beta = 1, nu = 1,
                      param = c(mu,delta,beta,nu),
                      main="Skew Hyperbolic Student-t QQ Plot",
                      xlab="Theoretical Quantiles",
                      ylab="Sample Quantiles",
                      plot.it = TRUE, line = TRUE,...){

    if (has.na <- any(ina <- is.na(y))) {
        yN <- y
        y <- y[!ina]
    }
    if (0 == (n <- length(y)))
        stop("y is empty or has only NAs")


    if (length(param) > 0){
        #check parameters
        parResult <- skewhypCheckPars(param)
        case <- parResult$case
        errMessage <- parResult$errMessage
        if(case == "error") stop(errMessage)
        mu <- param[1]
        delta <- param[2]
        beta <- param[3]
        nu <- param[4]
    } else {
        fitResults <- skewhypFit(y, freq = NULL, breaks = NULL,
            ParamStart = NULL, startMethod = "Nelder-Mead",
            startValues = "LA", method = "Nelder-Mead",
            hessian = FALSE, plots = FALSE, printOut = FALSE,
            controlBFGS = list(maxit = 200), controlNLM = list(maxit = 1000),
            maxitNLM = 1500)
        param <- fitResults$param
        names(param) = NULL
        mu <- param[1]
        delta <- param[2]
        beta <- param[3]
        nu <- param[4]
    }

    x <- qskewhyp(ppoints(n), param = param)[order(order(y))]

    if (has.na) {
        y <- x
        x <- yN
        x[!ina] <- y
        y <- yN
    }
    if (plot.it) {
        qqplot(x, y, main = main, xlab = xlab, ylab = ylab)
        title(sub = paste("param = (", round(param[1], 3), ",",
              round(param[2], 3), ",", round(param[3], 3), ",",
              round(param[4], 3), ")", sep = ""), ...)
    }
    if (line)
        abline(0, 1)
    invisible(list(x = x, y = y))

}


###### PP Plot function ######################################################
ppskewhyp <- function(y, beta = NULL, delta = NULL, mu = NULL, nu = NULL,
                      param = c(mu,delta,beta,nu),
                      main = "Skew Hyperbolic Student-t P-P Plot",
                      xlab = "Uniform Quantiles",
                      ylab = "Probability-integral-transformed Data",
                      plot.it = TRUE, line = TRUE, ...) {

    if (has.na <- any(ina <- is.na(y))) {
        yN <- y
        y <- y[!ina]
    }
    if (0 == (n <- length(y)))
        stop("data is empty")

    if (length(param) > 0) {
        #check parameters
        parResult <- skewhypCheckPars(param)
        case <- parResult$case
        errMessage <- parResult$errMessage
        if(case == "error") stop(errMessage)
        mu <- param[1]
        delta <- param[2]
        beta <- param[3]
        nu <- param[4]

    } else {
        fitResults <- skewhypFit(y, freq = NULL, breaks = NULL,
              paramStart = NULL, startMethod = "Nelder-Mead",
              startValues = "SL", method = "Nelder-Mead", hessian = FALSE,
              plots = FALSE, printOut = FALSE, controlBFGS = list(maxit = 200),
              controlNLM = list(maxit = 1000), maxitNLM = 1500, ...)

        param <- fitResults$param
        names(param) = NULL
        mu <- param[1]
        delta <- param[2]
        beta <- param[3]
        nu <- param[4]
    }

    yvals <- pskewhyp(y, param = param)
    xvals <- ppoints(n, a = 1/2)[order(order(y))]
    if (has.na) {
        y <- yvals
        x <- xvals
        yvals <- yN
        yvals[!ina] <- y
        xvals <- yN
        xvals[!ina] <- x
    }
    if (plot.it) {
        plot(xvals, yvals, main = main, xlab = xlab, ylab = ylab,
             ylim = c(0, 1), xlim = c(0, 1), ...)
        title(sub = paste("param = (", round(param[1], 3), ",",
              round(param[2], 3), ",", round(param[3], 3), ",",
              round(param[4], 3), ")", sep = ""), ...)
    }
    if (line)
        abline(0, 1)
    invisible(list(x = xvals, y = yvals))
}
