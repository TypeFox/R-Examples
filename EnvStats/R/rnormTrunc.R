rnormTrunc <-
function (n, mean = 0, sd = 1, min = -Inf, max = Inf) 
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
    arg.mat <- cbind.no.warn(n.vec = rep(1, n), mean = as.vector(mean), 
        sd = as.vector(sd), min = as.vector(min), max = as.vector(max))
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
        mean <- arg.mat[!na.index, "mean"]
        sd <- arg.mat[!na.index, "sd"]
        min <- arg.mat[!na.index, "min"]
        max <- arg.mat[!na.index, "max"]
        if (any(sd < .Machine$double.eps)) 
            stop("All non-missing values of 'sd' must be positive.")
        if (any(min > max)) 
            stop(paste("All non-missing values of 'min' must be", 
                "less than the corresponding elements of 'max'."))
        r[!na.index] <- qnormTrunc(p = runif(n.vec), mean = mean, 
            sd = sd, min = min, max = max)
        return(r)
    }
}
