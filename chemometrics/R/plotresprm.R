plotresprm <-
function (prmdcvobj, optcomp, y, X, ...)
{
    prm.cv <- prm_cv(X,y, a = optcomp, plot.opt=FALSE, ...)
    par(mfrow = c(1, 2))
    predcv <- prm.cv$predicted[, optcomp]
    preddcvall <- prmdcvobj$pred[,optcomp, ]
    preddcv <- apply(preddcvall, 1, mean)
    ylimits <- max(abs(preddcvall - drop(y)))
    ylimits <- sort(c(-ylimits, ylimits))
    plot(predcv, predcv - y, xlab = "Predicted y", ylab = "Residuals",
        cex.lab = 1.2, cex = 0.7, pch = 3, col = 1, ylim = ylimits, ...)
    title("Results from CV")
    abline(h = 0, lty = 1)
    plot(preddcv, preddcv - y, xlab = "Predicted y", ylab = "Residuals",
        cex.lab = 1.2, cex = 0.7, pch = 3, col = gray(0.6), type = "n",
        ylim = ylimits, ...)
    for (i in 1:ncol(preddcvall)) {
        points(preddcv, preddcvall[, i] - y, cex = 0.7, pch = 3,
            col = gray(0.6))
    }
    points(preddcv, preddcv - y, cex = 0.7, pch = 3, col = 1)
    title("Results from Repeated Double-CV")
    abline(h = 0, lty = 1)
    invisible()
}
