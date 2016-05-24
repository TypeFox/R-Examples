plnormMix <-
function (q, meanlog1 = 0, sdlog1 = 1, meanlog2 = 0, sdlog2 = 1, 
    p.mix = 0.5) 
{
    names.q <- names(q)
    arg.mat <- cbind.no.warn(q = as.vector(q), meanlog1 = as.vector(meanlog1), 
        sdlog1 = as.vector(sdlog1), meanlog2 = as.vector(meanlog2), 
        sdlog2 = as.vector(sdlog2), p.mix = as.vector(p.mix))
    for (i in c("q", "meanlog1", "sdlog1", "meanlog2", "sdlog2", 
        "p.mix")) assign(i, arg.mat[, i])
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        p <- rep(NA, length(q))
    else {
        if (any(c(sdlog1[!na.index], sdlog2[!na.index]) < .Machine$double.eps)) 
            stop("All values of 'sdlog1' and 'sdlog2' must be positive.")
        if (any(p.mix[!na.index] < 0 | p.mix[!na.index] > 1)) 
            stop("All values of 'p.mix' must be between 0 and 1.")
        p <- (1 - p.mix) * plnorm(q, meanlog1, sdlog1) + p.mix * 
            plnorm(q, meanlog2, sdlog2)
    }
    if (!is.null(names.q)) 
        names(p) <- rep(names.q, length = length(p))
    else names(p) <- NULL
    p
}
