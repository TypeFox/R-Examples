"stats2"<-function (x, by,digits=8){
    if (!missing(by)) {
        x <- cat.to.list(c(x), by)
    }
    if (!is.list(x) & !is.matrix(x))
        x <- matrix(x, ncol = 1)
    if (is.list(x)) {
        ncol <- length(names(x))
        out <- matrix(NA, ncol = ncol, nrow = length(describe2()))
        dimnames(out) <- list(describe2(), names(x))
        for (j in (1:ncol)) {
            if (is.numeric(x[[j]])) {
                out[, j] <- describe2(x[[j]])
            }
        }
        return(round(out,digits=digits))
    }
    if (is.matrix(x)) {
        nc <- ncol(x)
        ex.skew<-rep(NA,nc)
        ex.kurt<-rep(NA,nc)
        out <- matrix(NA, ncol = nc, nrow = length(describe2()))
        dimnames(out) <- list(describe2(), dimnames(x)[[2]])
        for (j in (1:nc)) {
            out[, j] <- describe2(x[, j])
        }
        return(round(out,digits=digits))
    }
}
