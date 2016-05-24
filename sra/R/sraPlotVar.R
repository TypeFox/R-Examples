sraPlotVar <-
function (srafit, series = levels(srafit$data$rep), legend = TRUE, 
    xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, pch = 1, 
    ...) 
{
    ymin <- 0
    ymax <- max(srafit$data$var)
    if (is.null(ylim)) {
        ylim <- c(0, ymax + (ymax - ymin)/5)
    }
    if (is.null(xlim)) {
        xlim <- c(0, max(srafit$data$gen))
    }
    if (is.null(xlab)) {
        xlab <- "Time (generations)"
    }
    if (is.null(ylab)) {
        ylab <- "Var(phenotype)"
    }
    plot(NULL, NULL, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, 
        axes = FALSE, yaxs = "i", ...)
    axis(1)
    axis(2)
    for (r in as.character(series)) {
        points(srafit$data$gen[srafit$data$rep == r], srafit$data$var[srafit$data$rep == 
            r], type = "b", lty = 3, pch = pch)
        lines(srafit$data$gen[srafit$data$rep == r], srafit$pred[[r]]$var)
        lines(srafit$data$gen[srafit$data$rep == r], srafit$pred[[r]]$varA, 
            lty = 2)
    }
    location <- "topleft"
    if (legend) {
        sraPlotlegend(labels = names(coef(srafit)), estimates = coef(srafit), 
            AIC = AIC(srafit), location = location)
    }
}
