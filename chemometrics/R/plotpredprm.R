plotpredprm <-
function (prmdcvobj, optcomp, y, X, ...)
{
    prm.cv = prm_cv(X,y, a = optcomp, plot.opt=FALSE, ...)
    par(mfrow = c(1, 2))
    ylimits = range(prmdcvobj$pred[, optcomp, ])
    plot(y, prm.cv$predicted[, optcomp], xlab = "Measured y",
        ylab = "Predicted y", cex.lab = 1.2, cex = 0.7, pch = 3,
        col = 1, ylim = ylimits, ...)
    title("Prediction from CV")
    abline(c(0, 1), lty = 1)
    plot(y, apply(prmdcvobj$pred[, optcomp, ], 1, mean), xlab = "Measured y",
        ylab = "Predicted y", cex.lab = 1.2, type = "n", ylim = ylimits, ...)
    for (i in 1:ncol(prmdcvobj$pred[, optcomp, ])) {
        points(y, prmdcvobj$pred[, optcomp,i], pch = 3, cex = 0.7,
            col = gray(0.6))
    }
    points(y, apply(prmdcvobj$pred[, optcomp, ], 1, mean),
        pch = 3, cex = 0.7, col = 1)
    title("Prediction from Repeated Double-CV")
    abline(c(0, 1), lty = 1)
    invisible()
}
