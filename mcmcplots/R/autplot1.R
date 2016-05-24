autplot1 <- function(x, chain=1, lag.max=NULL, partial=FALSE, col=mcmcplotsPalette(1), style=c("gray", "plain"), ylim=NULL, ...){
    style <- match.arg(style)
    if (partial){
        ylab <- "Partial Autocorrelation"
        xacf <-  pacf(as.ts(x[[chain]]), lag.max = lag.max, plot = FALSE)
    } else {
        ylab <- "Autocorrelation"
        xacf <-  acf(as.ts(x[[chain]]), lag.max = lag.max, plot = FALSE)
    }
    clim <- c(-1, 1)*qnorm(0.975)/sqrt(xacf$n.used)
    for (j in 1:nvar(x)) {
        if (is.null(ylim)){
            ylim <- range(c(clim, xacf$acf[, j, j]))
        }
        if (style=="gray"){
            plot(xacf$lag[, j, j], xacf$acf[, j, j], type = "n", ylab = ylab, xlab = "Lag", ylim = ylim, bty="n", xaxt="n", yaxt="n", ...)
            .graypr()
            rect(par("usr")[1], clim[1], par("usr")[2], clim[2], col=rgb(0.5, 0.5, 0.5, 0.35), border=NA)
            lines(xacf$lag[, j, j], xacf$acf[, j, j], type="h", lwd=2, col=col)
        }
        if (style=="plain"){
            plot(xacf$lag[, j, j], xacf$acf[, j, j], type = "h", ylab = ylab, xlab = "Lag", ylim = ylim, lwd=2, ...)
            abline(h=c(0, clim), col="gray")
        }
    }
}
