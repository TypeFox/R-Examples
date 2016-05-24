dpareto <-
function (x, location, shape = 1) 
{
    names.x <- names(x)
    arg.mat <- cbind.no.warn(x = as.vector(x), location = as.vector(location), 
        shape = as.vector(shape))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        y <- rep(NA, nrow(arg.mat))
    else {
        y <- numeric(nrow(arg.mat))
        y[na.index] <- NA
        y.no.na <- y[!na.index]
        for (i in c("x", "location", "shape")) assign(i, arg.mat[!na.index, 
            i])
        if (any(c(location, shape) < .Machine$double.eps)) 
            stop("All values of 'location' and 'shape' must be positive.")
        x.out <- x < location
        y.no.na[x.out] <- 0
        if (any(index <- !x.out)) {
            for (i in c("x", "location", "shape")) assign(i, 
                get(i)[index])
            y.no.na[index] <- shape * (location^shape) * (x^(-(shape + 
                1)))
        }
        y[!na.index] <- y.no.na
    }
    if (!is.null(names.x)) 
        names(y) <- rep(names.x, length = length(y))
    else names(y) <- NULL
    y
}
