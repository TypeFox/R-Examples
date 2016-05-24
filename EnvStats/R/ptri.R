ptri <-
function (q, min = 0, max = 1, mode = 1/2) 
{
    names.q <- names(q)
    arg.mat <- cbind.no.warn(q = as.vector(q), min = as.vector(min), 
        max = as.vector(max), mode = as.vector(mode))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        p <- rep(NA, nrow(arg.mat))
    else {
        p <- numeric(nrow(arg.mat))
        p[na.index] <- NA
        p.no.na <- p[!na.index]
        for (i in c("q", "min", "max", "mode")) assign(i, arg.mat[!na.index, 
            i])
        if (any(is.infinite(min)) || any(is.infinite(max))) 
            stop("All non-missing values of 'min' and 'max' must be finite.")
        if (any(mode <= min) || any(max <= mode)) 
            stop(paste("All values of 'mode' must be larger than", 
                "the corresponding values of 'min', and all", 
                "values of 'max' must be larger than the", "corresponding values of 'mode'."))
        q.low <- q <= min
        p.no.na[q.low] <- 0
        q.high <- q >= max
        p.no.na[q.high] <- 1
        if (any(index <- !(q.low | q.high))) {
            for (i in c("q", "min", "max", "mode")) assign(i, 
                get(i)[index])
            mmm <- max - min
            p.no.na[index] <- ifelse(q <= mode, (q - min)^2/(mmm * 
                (mode - min)), 1 - ((max - q)^2/(mmm * (max - 
                mode))))
        }
        p[!na.index] <- p.no.na
    }
    if (!is.null(names.q)) 
        names(p) <- rep(names.q, length = length(p))
    else names(p) <- NULL
    p
}
