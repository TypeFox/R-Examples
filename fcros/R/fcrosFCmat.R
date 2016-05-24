fcrosFCmat <-
function(xdata, cont, test, log2.opt = 0, trim.opt = 0.25) {
    n <- nrow(xdata);
    idnames <- xdata[,1]; # first column is unique ID
    xcol <- colnames(xdata);
    n.xcol <- length(xcol);
    idx1 <- xcol %in% cont;
    m1 <- sum(idx1 == TRUE);
    idx2 <- xcol %in% test;
    m2 <- sum(idx2 == TRUE);
    m <- m1+m2;
    m1m2 <- m1*m2;

    # form fcros data matrix
    dmat <- matrix(c(rep(0,n*m)), ncol = m);
    x1 <- matrix(c(rep(0,n*m1)), ncol = m1);
    x2 <- matrix(c(rep(0,n*m2)), ncol = m2);
    if (log2.opt) {
         x1 <- log2(xdata[,idx1 == TRUE]);
    } else {
         x1 <- xdata[,idx1 == TRUE];
    }
    dmat[,1:m1] <- as.matrix(x1);
    if (log2.opt) {
         x2 <- log2(xdata[,idx2 == TRUE]);
    } else {
         x2 <- xdata[,idx2 == TRUE];
    }
    dmat[,(m1+1):m] <- as.matrix(x2);

    # compute matrix containing pairwise fold changes
    fcMat <- matrix(c(rep(0,n*m1m2)), ncol = m1m2);
    k <- 1;
    for (i in 1:m1) {
        for (j in 1:m2) {
            fcMat[,k] <- dmat[,m1+j]-dmat[,i];
            k <- k+1;
        }
    }
    # compute the fold changes
    FC <- matrix(c(rep(0,n)), ncol = 1);
    for (i in 1:n) {
        x1 <- dmat[i,1:m1];
        x2 <- dmat[i,(m1+1):m];
        FC[i] <- mean(2^x2)/mean(2^x1);
    }
    FC2 = apply(2^fcMat, 1, mean, trim = trim.opt);

    list(idnames=idnames, fcMat = fcMat, FC = FC, FC2 = FC2);
}
