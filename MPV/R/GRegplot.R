GRegplot <- function(X, y, sortTrt = FALSE, includeIntercept = TRUE, type=c("hist", "dot")) 
    {
    Qy.out <- GFplot(X, y, plotIt = FALSE, sortTrt, type = "hist", includeIntercept, 
labels)
    p <- ncol(X)
    n <- nrow(X)
    if (!("x0" %in% names(X))) {
        X <- cbind(x0=rep(1,n), X)
        p <- p+1
    }
    y <- Qy.out$treatments[1:min(n,p)]
    x <- Qy.out$errors
    n <- length(y)
    x <- c(y,x)
    IQRange <- IQR(x)
    s <- IQRange/1.34898
    z <- x[abs(x) < 3.5*s]
    s <- sd(z)
    m <- length(z)
    B <- s*qt(.5+.95^(1/n)/2, df=m-1)
    hor.range <- range(c(abs(x),B))*1.075
    if (type=="hist") {
        hist(abs(y), xlab=paste("Estimates of R", expression(beta)), 
            main=" ", xlim=hor.range)
        rug(B, col=2, lwd=2)
    } else {
        stripchart(abs(y), xlim=hor.range, xlab=expression(Q^T*y))
        points(B, 1, col=2, pch=3, cex=2)      
    }
}
