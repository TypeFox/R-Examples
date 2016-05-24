sraPlotMean <-
function (srafit, series = levels(srafit$data$rep), legend = TRUE, 
    xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, pch = 1, 
    ...) 
{
    ymin <- min(srafit$data$mean)
    ymax <- max(srafit$data$mean)
    if (is.null(ylim)) {
        ylim <- c(ymin - (ymax - ymin)/5, ymax + (ymax - ymin)/5)
    }
    if (is.null(xlim)) {
        xlim <- c(0, max(srafit$data$gen))
    }
    if (is.null(xlab)) {
        xlab <- "Time (generations)"
    }
    if (is.null(ylab)) {
        ylab <- "Phenotype"
    }
    plot(NULL, NULL, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, 
        axes = FALSE, ...)
    axis(1)
    axis(2)
    for (r in as.character(series)) {
        points(srafit$data$gen[srafit$data$rep == r], srafit$data$mean[srafit$data$rep == 
            r], type = "b", lty = 3, pch = pch)
        lines(srafit$data$gen[srafit$data$rep == r], srafit$pred[[r]]$phen)
    }
    location <- if (ymax - srafit$pred[[1]]$phen[1] > (ymax - 
        ymin)/2) 
        "topleft"
    else "bottomleft"
    if (legend) {
        sraPlotlegend(labels = names(coef(srafit)), estimates = coef(srafit), 
            AIC = AIC(srafit), location = location)
    }
}
