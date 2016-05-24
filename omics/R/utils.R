`%din%` <- function(x, y) {
    by <- intersect(names(x), names(y))
    nx <- nrow(x <- as.data.frame(x))
    ny <- nrow(y <- as.data.frame(y))
    bx <- x[,by,drop=FALSE]
    by <- y[,by,drop=FALSE]
    names(bx) = names(by) <- paste("V", seq_len(ncol(bx)), sep="")
    bz <- do.call(paste, c(rbind(bx, by), sep="\r"))
    bx <- bz[seq_len(nx)]
    by <- bz[nx + seq_len(ny)]
    comm <- match(bx, by, 0)
    x[comm > 0,]
}

na.count <- function(X, margin, fraction=TRUE) {
    result <- apply(is.na(X), margin, sum)
    if (fraction) {
        result <- result / prod(dim(X)[-margin])
    }
    result
}

re.match <- function(pattern, x, ...) {
    do.call(rbind, regmatches(x, regexec(pattern, x, ...)))
}
