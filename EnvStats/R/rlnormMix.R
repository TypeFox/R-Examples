rlnormMix <-
function (n, meanlog1 = 0, sdlog1 = 1, meanlog2 = 0, sdlog2 = 1, 
    p.mix = 0.5) 
{
    ln <- length(n)
    if (ln == 0) 
        stop("'n' must be a non-empty scalar or vector.")
    if (ln > 1) 
        n <- ln
    else {
        if (is.na(n) || n <= 0 || n != trunc(n)) 
            stop("'n' must be a positive integer or a vector.")
    }
    if (length(p.mix) != 1 || is.na(p.mix) || p.mix < 0 || p.mix > 
        1) 
        stop("'p.mix' must be a single number between 0 and 1.")
    arg.mat <- cbind.no.warn(dum = rep(1, n), meanlog1 = as.vector(meanlog1), 
        sdlog1 = as.vector(sdlog1), meanlog2 = as.vector(meanlog2), 
        sdlog2 = as.vector(sdlog2))[, -1, drop = FALSE]
    if (n < nrow(arg.mat)) 
        arg.mat <- arg.mat[1:n, , drop = FALSE]
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        return(rep(NA, n))
    else {
        r <- numeric(n)
        r[na.index] <- NA
        r.no.na <- r[!na.index]
        for (i in c("meanlog1", "sdlog1", "meanlog2", "sdlog2")) assign(i, 
            arg.mat[!na.index, i])
        if (any(c(sdlog1, sdlog2) < .Machine$double.eps)) 
            stop("All non-missing values of 'sdlog1' and 'sdlog2' must be positive.")
        n.no.na <- sum(!na.index)
        index <- rbinom(n.no.na, 1, 1 - p.mix)
        n1 <- sum(index)
        n2 <- n.no.na - n1
        if (n1 > 0) 
            r.no.na[index == 1] <- rlnorm(n1, meanlog1, sdlog1)
        if (n2 > 0) 
            r.no.na[index == 0] <- rlnorm(n2, meanlog2, sdlog2)
        r[!na.index] <- r.no.na
        return(r)
    }
}
