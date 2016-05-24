acf <-
function (x, lag.max = NULL, type = c("correlation", "covariance", 
    "partial")[1], plot = TRUE, na.action = na.fail, demean = TRUE, 
    drop.lag.0 = TRUE,  ...) 
{
    acf.out <- stats:::acf(x=x, lag.max=lag.max, type=type, plot=F, na.action=na.action,
        demean=demean,...)
    acf.out$series <- deparse(substitute(x))
    if (drop.lag.0) {
        if (type == "correlation") {
            acf.out$acf = acf.out$acf[-1, , , drop = FALSE]
            acf.out$lag = acf.out$lag[-1, , , drop = FALSE]
        }
    }
        if (plot) {
        plot1.acf(acf.out, ...)
        return(invisible(acf.out))
    }
    else return(acf.out)
}
