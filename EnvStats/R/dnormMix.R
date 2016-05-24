dnormMix <-
function (x, mean1 = 0, sd1 = 1, mean2 = 0, sd2 = 1, p.mix = 0.5) 
{
    names.x <- names(x)
    arg.mat <- cbind.no.warn(x = as.vector(x), mean1 = as.vector(mean1), 
        sd1 = as.vector(sd1), mean2 = as.vector(mean2), sd2 = as.vector(sd2), 
        p.mix = as.vector(p.mix))
    for (i in c("x", "mean1", "sd1", "mean2", "sd2", "p.mix")) assign(i, 
        arg.mat[, i])
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        y <- rep(NA, length(x))
    else {
        if (any(c(sd1[!na.index], sd2[!na.index]) < .Machine$double.eps)) 
            stop("All values of 'sd1' and 'sd2' must be positive.")
        if (any(p.mix[!na.index] < 0 | p.mix[!na.index] > 1)) 
            stop("All values of 'p.mix' must be between 0 and 1.")
        y <- (1 - p.mix) * dnorm(x, mean1, sd1) + p.mix * dnorm(x, 
            mean2, sd2)
    }
    if (!is.null(names.x)) 
        names(y) <- rep(names.x, length = length(y))
    else names(y) <- NULL
    y
}
