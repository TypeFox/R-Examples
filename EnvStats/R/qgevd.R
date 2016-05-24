qgevd <-
function (p, location = 0, scale = 1, shape = 0) 
{
    names.p <- names(p)
    arg.mat <- cbind.no.warn(p = as.vector(p), location = as.vector(location), 
        scale = as.vector(scale), shape = as.vector(shape))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        q <- rep(NA, nrow(arg.mat))
    else {
        q <- numeric(nrow(arg.mat))
        q[na.index] <- NA
        q.no.na <- q[!na.index]
        for (i in c("p", "location", "scale", "shape")) assign(i, 
            arg.mat[!na.index, i])
        if (any(p < 0 | p > 1)) 
            stop("All values of 'p' must be between 0 and 1.")
        if (any(scale < .Machine$double.eps)) 
            stop("All values of 'scale' must be positive.")
        q.no.na[p == 0] <- -Inf
        q.no.na[p == 1] <- Inf
        index <- (0 < p) & (p < 1)
        ind <- index & shape == 0
        if (any(ind)) 
            q.no.na[ind] <- location[ind] - scale[ind] * log(-log(p[ind]))
        ind <- index & shape != 0
        if (any(ind)) 
            q.no.na[ind] <- location[ind] + (scale[ind] * (1 - 
                (-log(p[ind]))^shape[ind]))/shape[ind]
        q[!na.index] <- q.no.na
    }
    if (!is.null(names.p)) 
        names(q) <- rep(names.p, length = length(q))
    else names(q) <- NULL
    q
}
