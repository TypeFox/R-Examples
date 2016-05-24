qevd <-
function (p, location = 0, scale = 1) 
{
    names.p <- names(p)
    arg.mat <- cbind.no.warn(p = as.vector(p), location = as.vector(location), 
        scale = as.vector(scale))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        q <- rep(NA, length(p))
    else {
        q <- numeric(nrow(arg.mat))
        q[na.index] <- NA
        q.no.na <- q[!na.index]
        for (i in c("p", "location", "scale")) assign(i, arg.mat[!na.index, 
            i])
        if (any(p < 0 | p > 1)) 
            stop("All non-missing values of 'p' must be between 0 and 1.")
        if (any(scale < .Machine$double.eps)) 
            stop("All non-missing values of 'scale' must be positive.")
        q.no.na[p == 0] <- -Inf
        q.no.na[p == 1] <- Inf
        index <- (0 < p) & (p < 1)
        q.no.na[index] <- location[index] - scale[index] * log(log(1/p[index]))
        q[!na.index] <- q.no.na
    }
    if (!is.null(names.p)) 
        names(q) <- rep(names.p, length = length(q))
    q
}
