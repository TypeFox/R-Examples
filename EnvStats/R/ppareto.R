ppareto <-
function (q, location, shape = 1) 
{
    names.q <- names(q)
    arg.mat <- cbind.no.warn(q = as.vector(q), location = as.vector(location), 
        shape = as.vector(shape))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        p <- rep(NA, nrow(arg.mat))
    else {
        p <- numeric(nrow(arg.mat))
        p[na.index] <- NA
        p.no.na <- p[!na.index]
        for (i in c("q", "location", "shape")) assign(i, arg.mat[!na.index, 
            i])
        if (any(c(location, shape) < .Machine$double.eps)) 
            stop("All values of 'location' and 'shape' must be positive.")
        q.out <- q <= location
        p.no.na[q.out] <- 0
        if (any(index <- !q.out)) {
            for (i in c("q", "location", "shape")) assign(i, 
                get(i)[index])
            p.no.na[index] <- 1 - (location/q)^shape
        }
        p[!na.index] <- p.no.na
    }
    if (!is.null(names.q)) 
        names(p) <- rep(names.q, length = length(p))
    else names(p) <- NULL
    p
}
