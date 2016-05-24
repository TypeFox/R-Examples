#' Compute and Plot the Forecasts Based on a Fitted DDE Model
#'
#' The fitted process and forecasts are plotted.
#' These may be either plotted simultaneously, as \code{matplot} does for
#' multivariate data, or one by one with a mouse click to move from one
#' plot to another. The function also accepts the other \code{plot}
#' specification arguments that the regular plot does. The forecast
#' is done by integreting the process forward and forecast confidence
#' is based on an arima model fitted to the residual using
#' \code{auto.arima} function from \code{forecast} package.
#' @param y The observed data matrix
#' @param times A sequence of time points at which data are observed. must be evenly spaced to be able to fit arima model to the residual.
#' @param h How many lags to predict forward.
#' @param pars The fitted parameters for the DDE.
#' @param beta The fitted parameters for the contribution of lags of delays.
#' @param proc The \code{proc} object returned from estimation functions.
#' @param more An object specifying additional arguments to fn.
#' @param tau A list of delay lags.
#' @param ndelay A vector inidicating which state process has a delay term.
#' @param fdobj0 A functional data object fitted with the initial history part of the data.
#' @param fdobj.d A functional data object fitted by generalized profiling.
#' @param ask a logical value: If TRUE, each curve is shown separately, and the plot advances with a mouse click.
#' @param xlab a label for the horizontal axis.
#' @param ylab a label for the vertical axis.
#' @param xlim a vector of length 2 containing axis limits for the horizontal axis.
#' @param ylim a vector of length 2 containing axis limits for the vertical axis.
#' @param axes Either a logical or a list or \code{NULL}.
#' ##' \describe{
##' \item{logical}{whether axes should be drawn on the plot}
##' \item{list}{a list used to create custom axes used to create axes via     \code{do.call(x$axes[[1]], x$axes[-1])}.}
##' }
#' @param ... additional plotting arguments that can be used with function \code{plot}.
#'
#' @return 'done'
#' @export
#'@seealso \code{\link{IntegrateForward.DDE}}
forecast.DDE <- function(y, times, h, pars, beta, proc, more, tau, ndelay, fdobj0, fdobj.d, ask = FALSE, xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, axes = NULL,...){
    if (is.null(axes)) {
        if (is.null(fdobj.d$basis$axes)) {
            Axes <- TRUE
            axFun <- FALSE
        }
        else {
            if (!inherits(fdobj.d$basis$axes, "list"))
                stop("fdobj.d$basis$axes must be a list;  ",
                     "class(fdobj.d$basis$axes) = ", class(fdobj.d$basis$axes))
            if (!(inherits(fdobj.d$basis$axes[[1]], "character") ||
                  inherits(fdobj.d$basis$axes[[1]], "function")))
                stop("fdobj.d$basis$axes[[1]] must be either a function or the ",
                     "name of a function;  class(fdobj.d$basis$axes[[1]]) = ",
                     class(fdobj.d$basis$axes[[1]]))
            Axes <- FALSE
            axFun <- TRUE
            axList <- c(fdobj.d$basis$axes, ...)
        }
    }
    else {
        if (is.logical(axes)) {
            Axes <- axes
            axFun <- FALSE
        }
        else {
            if (!inherits(axes, "list"))
                stop("axes must be a logical or a list;  class(axes) = ",
                     class(axes))
            if (!(inherits(axes[[1]], "character") || inherits(axes[[1]],
                                                               "function")))
                stop("axes[[1]] must be either a function or the ",
                     "name of a function;  class(axes[[1]]) = ",
                     class(axes[[1]]))
            Axes <- FALSE
            axFun <- TRUE
            axList <- c(axes, ...)
        }
    }
    fdmat <- eval.fd(times, fdobj.d)
    basisobj <- fdobj.d$basis
    rangex <- range(times)
    times.fore <- rangex[2] + (0:h) * (rangex[2] - rangex[1]) / length(times)
    xlim <- c(rangex[1], max(times.fore))
    res <- y - fdmat
    forward.obj <- IntegrateForward.DDE( times.forecast = times.fore, pars = pars, beta = beta,  proc = proc, more = more, tau = tau, ndelay = ndelay, fdobj0 = fdobj0, fdobj.d = fdobj.d)$states
    if(!is.matrix(forward.obj)){
        forward.obj <- matrix(forward.obj, ncol = 1)
    }
    res.forecast <- list()
    for(i in 1:ncol(res)){
        points(times, y[,i])
        res.arima <- forecast::auto.arima(res[,i])
        res.forecast[[i]] <- forecast::forecast(res.arima, level = c(95), h = h)
        forward.obj[-1,i] <- res.forecast[[i]]$mean + forward.obj[-1,i]
    }
    fdmat <-rbind(fdmat, forward.obj[-1,])
    rangey <- range(fdmat)
       if (is.null(ylim))
            ylim <- rangey
    fdnames = fdobj.d$fdnames
    coef   <- fdobj.d$coefs
    coefd  <- dim(coef)
    nbasis    <- coefd[1]
    # Number of functional observations
    nrep   <- coefd[2]
    fdlabelslist = fdlabels(fdnames, nrep, 1)
    xlabel = fdlabelslist$xlabel
    ylabel = fdlabelslist$ylabel
    casenames = fdlabelslist$casenames
    varnames = fdlabelslist$varnames
    if (is.null(xlab))
        xlab <- xlabel
    if (is.null(ylab))
        ylab <- ylabel
    if(!ask){
        matplot(c(times, times.fore[-1]), fdmat, type = "l",  ask = FALSE, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, axes = Axes, ...)
        for(i in 1:ncol(res)){
            points(times, y[,i])
            upi <- res.forecast[[i]]$upper + forward.obj[-1,i] - res.forecast[[i]]$mean
            lpi <- res.forecast[[i]]$lower + forward.obj[-1,i] - res.forecast[[i]]$mean
            lines(times.fore[-1], upi, lty = 3)
            lines(times.fore[-1], lpi, lty = 3)
        }
    }
    else {
        op <- par(ask = FALSE)
        on.exit(par(op))
        cat("Multiple plots:  Click in the plot to advance to the next")
        for(i in 1:ncol(res)){
            plot(c(times, times.fore[-1]), fdmat[, i], type = "l", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, axes = Axes, ...)
            points(times, y[,i])
            if (!is.null(casenames)) title(casenames[i])
            else                     title(paste("Case",i))
            upi <- res.forecast[[i]]$upper + forward.obj[-1,i] - res.forecast[[i]]$mean
            lpi <- res.forecast[[i]]$lower + forward.obj[-1,i] - res.forecast[[i]]$mean
            lines(times.fore[-1], upi, lty = 3)
            lines(times.fore[-1], lpi, lty = 3)
        }
    }
    "done"
}

