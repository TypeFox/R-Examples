rnormMix <-
function (n, mean1 = 0, sd1 = 1, mean2 = 0, sd2 = 1, p.mix = 0.5) 
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
    arg.mat <- cbind.no.warn(dum = rep(1, n), mean1 = as.vector(mean1), 
        sd1 = as.vector(sd1), mean2 = as.vector(mean2), sd2 = as.vector(sd2))[, 
        -1, drop = FALSE]
    if (n < nrow(arg.mat)) 
        arg.mat <- arg.mat[1:n, , drop = FALSE]
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        return(rep(NA, n))
    else {
        r <- numeric(n)
        r[na.index] <- NA
        r.no.na <- r[!na.index]
        for (i in c("mean1", "sd1", "mean2", "sd2")) assign(i, 
            arg.mat[!na.index, i])
        if (any(c(sd1, sd2) < .Machine$double.eps)) 
            stop("All non-missing values of 'sd1' and 'sd2' must be positive.")
        n.no.na <- sum(!na.index)
        index <- rbinom(n.no.na, 1, 1 - p.mix)
        n1 <- sum(index)
        n2 <- n.no.na - n1
        if (n1 > 0) 
            r.no.na[index == 1] <- rnorm(n1, mean1, sd1)
        if (n2 > 0) 
            r.no.na[index == 0] <- rnorm(n2, mean2, sd2)
        r[!na.index] <- r.no.na
        return(r)
    }
}
