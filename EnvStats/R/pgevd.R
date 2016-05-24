pgevd <-
function (q, location = 0, scale = 1, shape = 0) 
{
    names.q <- names(q)
    arg.mat <- cbind.no.warn(q = as.vector(q), location = as.vector(location), 
        scale = as.vector(scale), shape = as.vector(shape))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        p <- rep(NA, nrow(arg.mat))
    else {
        p <- numeric(nrow(arg.mat))
        p[na.index] <- NA
        p.no.na <- p[!na.index]
        for (i in c("q", "location", "scale", "shape")) assign(i, 
            arg.mat[!na.index, i])
        if (any(scale < .Machine$double.eps)) 
            stop("All values of 'scale' must be positive.")
        z <- (q - location)/scale
        shape.eq.0 <- shape == 0
        if (any(shape.eq.0)) 
            p.no.na[shape.eq.0] <- exp(-exp(-z[shape.eq.0]))
        shape.gt.0 <- shape > 0
        if (any(shape.gt.0)) {
            z.out <- z > 1/shape
            p.no.na[shape.gt.0 & z.out] <- 1
            index <- shape.gt.0 & !z.out
            if (any(index)) 
                p.no.na[index] <- exp(-(1 - shape[index] * z[index])^(1/shape[index]))
        }
        shape.lt.0 <- shape < 0
        if (any(shape.lt.0)) {
            z.out <- z < 1/shape
            p.no.na[shape.lt.0 & z.out] <- 0
            index <- shape.lt.0 & !z.out
            if (any(index)) 
                p.no.na[index] <- exp(-(1 - shape[index] * z[index])^(1/shape[index]))
        }
        p[!na.index] <- p.no.na
    }
    if (!is.null(names.q)) 
        names(p) <- rep(names.q, length = length(p))
    else names(p) <- NULL
    p
}
