qpareto <-
function (p, location, shape = 1) 
{
    names.p <- names(p)
    arg.mat <- cbind.no.warn(p = as.vector(p), location = as.vector(location), 
        shape = as.vector(shape))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        q <- rep(NA, nrow(arg.mat))
    else {
        q <- numeric(nrow(arg.mat))
        q[na.index] <- NA
        q.no.na <- q[!na.index]
        for (i in c("p", "location", "shape")) assign(i, arg.mat[!na.index, 
            i])
        if (any(p < 0) || any(p > 1)) 
            stop("All non-missing values of 'p' must be between 0 and 1.")
        if (any(c(location, shape) < .Machine$double.eps)) 
            stop("All values of 'location' and 'shape' must be positive.")
        q.no.na[p == 0] <- location[p == 0]
        q.no.na[p == 1] <- Inf
        if (any(index <- 0 < p & p < 1)) {
            for (i in c("p", "location", "shape")) assign(i, 
                get(i)[index])
            q.no.na[index] <- location * (1 - p)^(-(1/shape))
        }
        q[!na.index] <- q.no.na
    }
    if (!is.null(names.p)) 
        names(q) <- rep(names.p, length = length(q))
    else names(q) <- NULL
    q
}
