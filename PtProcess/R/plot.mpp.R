plot.mpp <- function(x, log = FALSE, ...){
    ls... <- list(...)
    times <- seq(x$TT[1], x$TT[2], length.out=1000)
    times <- sort(c(times, x$data$time))
    params <- x$params
    gparams <- eval(x$gmap)
    y <- x$gif(x$data, times, params=gparams)
    if (log==TRUE) y <- log(y)
    if (is.null(ls...$ylab)){
        if (log==TRUE)
            ylab <- expression(paste("log ", lambda[g](t)))
        else
            ylab <- expression(lambda[g](t))
    }
    else ylab <- ls...$ylab
    if (is.null(ls...$xlab)) xlab <- expression(t)
    else xlab <- ls...$xlab
    plot(times, y, type="l", ylab=ylab, xlab=xlab,
         xlim=x$TT)
}

