Uplot <- function (X.qr, Xcolumn=1, ...) {
    if (class(X.qr) !="qr") {
        X <- X.qr
        p <- ncol(X)
        n <- nrow(X)
        if (!("x0" %in% names(X))) {
            X <- cbind(x0=rep(1,n), X)
            p <- p+1
        } 
        X.qr <- qr(X)
    }
    R <- qr.R(X.qr)
    Xcols <- Xcolumn
    if (is.character(Xcols)) Xcols <- which(colnames(R) %in% Xcolumn)
    m <- length(Xcolumn) 
    if (m > 5) {
        layout <- c(5,ceiling(m/5)) 
    } else {
        layout <- c(m, 1)
    }
    par(mfcol=layout)
    for (i in Xcols) {
        barplot(abs(R[i,]), axes=FALSE, xlab="", ylab=colnames(R)[i], ...)
    }
}
