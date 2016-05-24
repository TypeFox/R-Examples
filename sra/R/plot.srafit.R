plot.srafit <-
function (x, series = levels(x$data$rep), resid = FALSE, variance = FALSE, 
    ...) 
{
    srafit <- x
    if (resid) {
        layout(cbind(c(2, 1)), heights = c(0.3, 0.6))
    }
    par(mar = c(5, 4, 1, 2))
    if (variance) 
        sraPlotVar(srafit, series, ...)
    else sraPlotMean(srafit, series, ...)
    if (resid) {
        par(mar = c(1, 4, 4, 2))
        if (variance) 
            sraPlotVarResid(srafit, series)
        else sraPlotMeanResid(srafit, series)
    }
}
