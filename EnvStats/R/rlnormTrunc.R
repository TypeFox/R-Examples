rlnormTrunc <-
function (n, meanlog = 0, sdlog = 1, min = 0, max = Inf) 
{
    ln <- length(n)
    if (ln < 1) 
        stop("'n' must be non-empty.")
    if (ln > 1) 
        n <- ln
    else {
        if (is.na(n) || n <= 0 || n != trunc(n)) 
            stop("'n' must be a positive integer or vector.")
    }
    arg.mat <- cbind.no.warn(n.vec = rep(1, n), meanlog = as.vector(meanlog), 
        sdlog = as.vector(sdlog), min = as.vector(min), max = as.vector(max))
    if (n < nrow(arg.mat)) 
        arg.mat <- arg.mat[1:n, , drop = FALSE]
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        return(rep(NA, n))
    else {
        r <- numeric(n)
        r[na.index] <- NA
        r.no.na <- r[!na.index]
        n.vec <- arg.mat[!na.index, "n.vec"]
        meanlog <- arg.mat[!na.index, "meanlog"]
        sdlog <- arg.mat[!na.index, "sdlog"]
        min <- arg.mat[!na.index, "min"]
        max <- arg.mat[!na.index, "max"]
        if (any(sdlog < .Machine$double.eps)) 
            stop("All non-missing values of 'sdlog' must be positive.")
        if (any(min > max)) 
            stop(paste("All non-missing values of 'min' must be", 
                "less than the corresponding elements of 'max'."))
        if (any(min < 0)) 
            stop(paste("All non-missing values of 'min' must be", 
                "greater than or equal to 0."))
        r[!na.index] <- qlnormTrunc(p = runif(n.vec), meanlog = meanlog, 
            sdlog = sdlog, min = min, max = max)
        return(r)
    }
}
