## plot a principle curve in PCA space
plot.prcurve <- function(x, axes = 1:2,
                         scaling = 0,
                         segments = TRUE,
                         col = "red",
                         col.seg = "forestgreen",
                         lwd = 2, lwd.seg = 1, ...) {
    pred <- predict(x$ordination, x$s, type = "wa",
                    scaling = scaling)[,axes]
    scrs <- scores(x$ordination, display = "sites", scaling = scaling,
                   choices = axes)
    xlim <- range(scrs[,1], pred[,1])
    ylim <- range(scrs[,2], pred[,2])
    plot(x$ordination, display = "sites", scaling = scaling, type = "n",
         xlim = xlim, ylim = ylim, choices = axes, ...)
    points(scrs, ...)
    if(segments)
        segments(scrs[,1], scrs[,2], pred[,1], pred[,2],
                 col = col.seg, lwd = lwd.seg)
    lines(pred[x$tag, 1:2], lwd = lwd, col = col, ...)
}
