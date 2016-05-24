qtri <-
function (p, min = 0, max = 1, mode = 1/2) 
{
    names.p <- names(p)
    arg.mat <- cbind.no.warn(p = as.vector(p), min = as.vector(min), 
        max = as.vector(max), mode = as.vector(mode))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        q <- rep(NA, nrow(arg.mat))
    else {
        q <- numeric(nrow(arg.mat))
        q[na.index] <- NA
        q.no.na <- q[!na.index]
        for (i in c("p", "min", "max", "mode")) assign(i, arg.mat[!na.index, 
            i])
        if (any(p < 0) || any(p > 1)) 
            stop("All non-missing values of 'p' must be between 0 and 1.")
        if (any(is.infinite(min)) || any(is.infinite(max))) 
            stop("All non-missing values of 'min' and 'max' must be finite.")
        if (any(mode <= min) || any(max <= mode)) 
            stop(paste("All values of 'mode' must be larger than", 
                "the corresponding values of 'min', and all", 
                "values of 'max' must be larger than the", "corresponding values of 'mode'."))
        q.no.na[p == 0] <- min[p == 0]
        q.no.na[p == 1] <- max[p == 1]
        if (any(index <- 0 < p & p < 1)) {
            for (i in c("p", "min", "max", "mode")) assign(i, 
                get(i)[index])
            mmm <- max - min
            q.no.na[index] <- ifelse(p <= ptri(mode, min = min, 
                max = max, mode = mode), min + sqrt(mmm * (mode - 
                min) * p), max - sqrt(mmm * (max - mode) * (1 - 
                p)))
        }
        q[!na.index] <- q.no.na
    }
    if (!is.null(names.p)) 
        names(q) <- rep(names.p, length = length(q))
    else names(q) <- NULL
    q
}
