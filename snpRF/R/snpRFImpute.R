snpRFImpute <- function(x.autosome=NULL,x.xchrom=NULL,x.covar, y, iter=5, ntree=300, ...) {
    if (any(is.na(y))) stop("Can't have NAs in", deparse(substitute(y)))
    if (!is.null(x.autosome)) {
       if (any(is.na(x.autosome))) stop("Can't have NAs in", deparse(substitute(x.autosome)))
    }
    if (!is.null(x.xchrom)) {
       if (any(is.na(x.xchrom))) stop("Can't have NAs in", deparse(substitute(x.xchrom)))
    }
    if (!any(is.na(x.covar))) stop("No NAs found in ", deparse(substitute(x)))
    xf <- na.roughfix(x.covar)
    hasNA <- which(apply(x.covar, 2, function(x) any(is.na(x))))
    if (is.data.frame(x.covar)) {
        isfac <- sapply(x.covar, is.factor)
    } else {
        isfac <- rep(FALSE, ncol(x.covar))
    }
    
    for (i in 1:iter) {

        prox <- snpRF(x.autosome=x.autosome,x.xchrom=x.xchrom,x.covar=xf, y=y, ntree=ntree, ..., 
	       	      do.trace=ntree,proximity=TRUE)$proximity

        for (j in hasNA) {

            miss <- which(is.na(x.covar[, j]))
            if (isfac[j]) {
                lvl <- levels(x.covar[[j]])
                catprox <- apply(prox[-miss, miss, drop=FALSE], 2,
                                 function(v) lvl[which.max(tapply(v, x.covar[[j]][-miss], mean))])
                xf[miss, j] <- catprox
            } else {
                sumprox <- colSums(prox[-miss, miss, drop=FALSE])
                xf[miss, j] <- (prox[miss, -miss, drop=FALSE] %*% xf[,j][-miss]) / (1e-8 + sumprox)
            }
            NULL
        }
    }
    xf <- cbind(y, xf)
    names(xf)[1] <- deparse(substitute(y))
    xf
}
