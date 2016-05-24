"plt.var" =
function (res, d1, d2, int = 0.5, var.label = FALSE, 
                     Xnames = NULL, Ynames = NULL) 
{
    if (!var.label) {
        plot(0, type = "n", xlim = c(-1, 1), ylim = c(-1, 1), 
        xlab = paste("Dimension ", d1), ylab = paste("Dimension ", d2))
        points(res$scores$corr.X.xscores[, d1], 
               res$scores$corr.X.xscores[, d2], 
               pch = 20, cex = 1.2, col = "red")
        points(res$scores$corr.Y.xscores[, d1], 
               res$scores$corr.Y.xscores[, d2], 
               pch = 24, cex = 0.7, col = "blue")
    }
    else {
        if (is.null(Xnames)) Xnames = res$names$Xnames
        if (is.null(Ynames)) Ynames = res$names$Ynames
        plot(0, type = "n", xlim = c(-1, 1), ylim = c(-1, 1), 
        xlab = paste("Dimension ", d1), ylab = paste("Dimension ", d2))
        text(res$scores$corr.X.xscores[, d1], 
             res$scores$corr.X.xscores[, d2], 
             Xnames, col = "red", font = 2)
        text(res$scores$corr.Y.xscores[, d1], 
             res$scores$corr.Y.xscores[, d2], 
             Ynames, col = "blue", font = 3)
    }

    abline(v = 0, h = 0)
    lines(cos(seq(0, 2 * pi, l = 100)), sin(seq(0, 2 * pi, l = 100)))
    lines(int * cos(seq(0, 2 * pi, l = 100)), 
          int * sin(seq(0, 2 * pi, l = 100)))
}

