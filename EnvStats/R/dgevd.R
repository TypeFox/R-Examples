dgevd <-
function (x, location = 0, scale = 1, shape = 0) 
{
    names.x <- names(x)
    arg.mat <- cbind.no.warn(x = as.vector(x), location = as.vector(location), 
        scale = as.vector(scale), shape = as.vector(shape))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        y <- rep(NA, nrow(arg.mat))
    else {
        y <- numeric(nrow(arg.mat))
        y[na.index] <- NA
        y.no.na <- y[!na.index]
        for (i in c("x", "location", "scale", "shape")) assign(i, 
            arg.mat[!na.index, i])
        if (any(scale < .Machine$double.eps)) 
            stop("All values of 'scale' must be positive.")
        z <- (x - location)/scale
        shape.eq.0 <- shape == 0
        if (any(shape.eq.0)) 
            y.no.na[shape.eq.0] <- (exp(-exp(-z[shape.eq.0])) * 
                exp(-z[shape.eq.0]))/scale[shape.eq.0]
        shape.gt.0 <- shape > 0
        if (any(shape.gt.0)) {
            z.out <- z >= 1/shape
            y.no.na[shape.gt.0 & z.out] <- 0
            index <- shape.gt.0 & !z.out
            if (any(index)) 
                y.no.na[index] <- exp(-(1 - shape[index] * z[index])^(1/shape[index])) * 
                  (1/scale[index]) * ((1 - shape[index] * z[index])^((1/shape[index]) - 
                  1))
        }
        shape.lt.0 <- shape < 0
        if (any(shape.lt.0)) {
            z.out <- z <= 1/shape
            y.no.na[shape.lt.0 & z.out] <- 0
            index <- shape.lt.0 & !z.out
            if (any(index)) 
                y.no.na[index] <- exp(-(1 - shape[index] * z[index])^(1/shape[index])) * 
                  (1/scale[index]) * ((1 - shape[index] * z[index])^((1/shape[index]) - 
                  1))
        }
        y[!na.index] <- y.no.na
    }
    if (!is.null(names.x)) 
        names(y) <- rep(names.x, length = length(x))
    else names(y) <- NULL
    y
}
