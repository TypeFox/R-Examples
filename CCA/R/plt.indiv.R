"plt.indiv" <-
function (res, d1, d2, ind.names = NULL) 
{
    if (is.null(ind.names)) ind.names = res$names$ind.names
    if (is.null(ind.names)) ind.names = 1:nrow(res$scores$xscores)
    plot(res$scores$xscores[, d1], res$scores$xscores[, d2], 
        type = "n", main = "", xlab = paste("Dimension ", d1), 
        ylab = paste("Dimension ", d2))
    text(res$scores$xscores[, d1], res$scores$xscores[, d2], 
        ind.names)
    abline(v = 0, h = 0, lty=2)
}

