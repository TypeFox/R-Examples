plotSEPprm <- function (prmdcvobj, optcomp, y, X, complete=TRUE, ...)
{
    matplot(prmdcvobj$SEPtrim, type = "l", lty = 1, col = gray(0.6),
        xlab = "Number of components", ylab = "Trimmed SEP", cex.lab = 1.2, ...)
    if (complete) {
        ncomp <- nrow(prmdcvobj$SEPtrim)
    }
    else {
        ncomp <- optcomp
    }
    prm.cv = prm_cv(X,y, a = ncomp, plot.opt=FALSE, ...)
    lines(prm.cv$SEPtrim, col = 1, lwd = 2)
    abline(v = optcomp, lty = 2)
    abline(h = prmdcvobj$SEPopt, lty = 2)
    list(SEPdcv = t(prmdcvobj$SEPtrim), SEPcv = prm.cv$SEPtrim)
}
