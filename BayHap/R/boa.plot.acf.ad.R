`boa.plot.acf.ad` <-
function (x, pname, annotate = boa.par("legend")) 
{
    boa.init()
    chain.import(x)
    boa.par(acf.lags=1)

   
    drawn <- FALSE
    parm <- boa:::boa.getparms(x, pname)
    if (is.matrix(parm)) {
        drawn <- TRUE
        val <- boa.par("par")
        cex <- ifelse(is.null(val$cex), 1, val$cex)
        lwd <- ifelse(is.null(val$lwd), 1, val$lwd)
        result <- acf(parm, plot = FALSE)
        plot(result$lag, result$acf, xlab = "Lag", ylab = "Autocorrelation", 
            main = pname, ylim = c(-1, 1), type = "h", lwd = lwd)
        abline(0, 0)
        usr <- par("usr")
        if (annotate) 
            legend(x = usr[2], y = 1, xjust = 1, yjust = 1, cex = cex, 
                legend = substring("x", first = 1, last = 16), 
                bty = "n")
    }
    return(drawn)
}

