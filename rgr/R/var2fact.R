var2fact <-
function (xmat) 
{
    ncols <- ncol(xmat)
    nrows <- nrow(xmat)
    names <- dimnames(xmat)[[2]]
    nind <- nrows * ncols
    xx.1 <- character(nind)
    xx.2 <- numeric(nind)
    for (i in 1:ncols) {
        i1 <- (i - 1) * nrows + 1
        i2 <- i1 + nrows - 1
        xx.1[i1:i2] <- names[i]
        xx.2[i1:i2] <- xmat[1:nrows, i]
    }
    xx <- cbind(xx.1, xx.2)
    dimnames(xx)[[2]] <- c("xx.1", "xx.2")
    invisible(xx)
}
