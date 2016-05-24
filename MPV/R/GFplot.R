GFplot <- function (X, y, plotIt=TRUE, sortTrt=FALSE, type="hist", includeIntercept=TRUE, labels=FALSE) {
    p <- ncol(X)
    n <- nrow(X)
    if (!("x0" %in% names(X))) {
        X <- cbind(x0=rep(1,n), X)
        p <- p+1
    } 
    X.qr <- qr(X)
    Q <- qr.Q(X.qr, complete=TRUE)
    Qty <- t(Q)%*%y
    Rb <- Qty[1:p]
    if (sortTrt) {
        newOrder <- rev(order(abs(Rb)))
        X <- X[, newOrder]
        outSortF <- GFplot(X, y, plotIt=FALSE, sortTrt=FALSE)
        errors <- outSortF$errors
        Rb <- outSortF$treatments
    } else {
        errors <- Qty[-(1:p)]
    }
    if (plotIt) {
        if (!includeIntercept) {
            index <- which(names(X)=="x0")
            Rb <- Rb[-index]
            X <- X[,-index]
            p <- p-1
        } 
        if (type=="QQ") {
            plotregionmin <- min(Rb, errors)
            plotregionmax <- max(Rb, errors)
            qqANOVA(errors, Rb, xlim=range(c(Rb, errors)),
               ylim=range(c(Rb, errors)), ylab="predictor effects")
            abline(0,1)
            if (p > 1) {            
                legend(plotregionmin, plotregionmax, legend=paste(names(X), 
                       c(paste("+", names(X)[-c(1,p)], "+ ..."), paste("+", names(X)[p] ),  " ")), 
                       col=as.numeric(factor(Rb)), pch=16)
            }
        } else {
            hor.range <- range(abs(Rb), abs(errors))*1.1
            hist(abs(errors), xlab="errors", main=" ", xlim=hor.range)
            rug(abs(Rb), col="red", lwd=3)
            if (labels) {
                 mtext(names(X), at=abs(Rb), side=1)
            } else {
                 if (sortTrt) {
                     mtext(newOrder+!includeIntercept, at=abs(Rb), side=1)
                 } 
                 else
                 {
                     mtext(1:ncol(X)+!includeIntercept, at=abs(Rb), side=1)
                 }
            }
        }
    } else {
        regF <- list(treatments=Rb, errors=errors, X=X, y=y, QR=X.qr)
        regF
    }
}
