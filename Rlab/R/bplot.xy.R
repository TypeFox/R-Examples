"bplot.xy" <-
function (x, y, N = 10, breaks = pretty(x, N,eps.correct=1), style = "tukey", 
    outlier = TRUE, plot = TRUE, xaxt = "s", ...) 
{
    out <- list()
    NBIN <- length(breaks) - 1
    centers <- (breaks[1:NBIN] + breaks[2:(NBIN + 1)])/2
    obj <- as.list(1:NBIN)
    names(obj) <- format(1:NBIN)
    for (k in 1:NBIN) {
        obj[[k]] <- describe.bplot(y[x < breaks[k + 1] & x > 
            breaks[k]], style = style, outlier = outlier)
    }
    if (plot) {
        bplot.obj(obj, pos = centers, label.cex = 0, outlier = outlier, 
            , xaxt = xaxt, ...)
    }
    else {
        return(list(centers = centers, breaks = breaks, bplot.obj = obj))
    }
    invisible()
}
