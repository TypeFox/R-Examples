dtri <-
function (x, min = 0, max = 1, mode = 1/2) 
{
    names.x <- names(x)
    arg.mat <- cbind.no.warn(x = as.vector(x), min = as.vector(min), 
        max = as.vector(max), mode = as.vector(mode))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        y <- rep(NA, nrow(arg.mat))
    else {
        y <- numeric(nrow(arg.mat))
        y[na.index] <- NA
        y.no.na <- y[!na.index]
        for (i in c("x", "min", "max", "mode")) assign(i, arg.mat[!na.index, 
            i])
        if (any(is.infinite(min)) || any(is.infinite(max))) 
            stop("All non-missing values of 'min' and 'max' must be finite.")
        if (any(mode <= min) || any(max <= mode)) 
            stop(paste("All values of 'mode' must be larger than", 
                "the corresponding values of 'min', and all", 
                "values of 'max' must be larger than the", "corresponding values of 'mode'."))
        mmm <- max - min
        y.no.na <- 2 * ifelse(x <= mode, (x - min)/(mmm * (mode - 
            min)), (max - x)/(mmm * (max - mode)))
        y.no.na[y.no.na < 0] <- 0
        y[!na.index] <- y.no.na
    }
    if (!is.null(names.x)) 
        names(y) <- rep(names.x, length = length(y))
    else names(y) <- NULL
    y
}
