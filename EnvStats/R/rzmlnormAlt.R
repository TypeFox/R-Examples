rzmlnormAlt <-
function (n, mean = exp(1/2), cv = sqrt(exp(1) - 1), p.zero = 0.5) 
{
    ln <- length(n)
    if (ln < 1) 
        stop("'n' must be non-empty.")
    if (ln > 1) 
        n <- ln
    else {
        if (is.na(n) || n <= 0 || n != trunc(n)) 
            stop("'n' must be a positive integer or a vector.")
    }
    if (length(p.zero) != 1 || is.na(p.zero) || p.zero <= 0 || 
        p.zero >= 1) 
        stop("'p.zero' must be a single number greater than 0 and less than 1.")
    arg.mat <- cbind.no.warn(dum = rep(1, n), mean = as.vector(mean), 
        cv = as.vector(cv))[, -1, drop = FALSE]
    if (n < nrow(arg.mat)) 
        arg.mat <- arg.mat[1:n, , drop = FALSE]
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        return(rep(NA, n))
    else {
        r <- numeric(n)
        r[na.index] <- NA
        r.no.na <- r[!na.index]
        for (i in c("mean", "cv")) assign(i, arg.mat[!na.index, 
            i])
        if (any(c(mean, cv) < .Machine$double.eps)) 
            stop("All non-missing values of 'mean' and 'cv' must be positive.")
        n.no.na <- sum(!na.index)
        index <- rbinom(n.no.na, 1, p.zero)
        n1 <- sum(index)
        n2 <- n.no.na - n1
        r.no.na[index == 1] <- 0
        if (n2 > 0) 
            r.no.na[index == 0] <- rlnormAlt(n2, mean, cv)
        r[!na.index] <- r.no.na
        return(r)
    }
}
