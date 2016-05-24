plnormAlt <-
function (q, mean = exp(1/2), cv = sqrt(exp(1) - 1), lower.tail = TRUE, 
    log.p = FALSE) 
{
    names.q <- names(q)
    arg.mat <- cbind.no.warn(q = as.vector(q), mean = as.vector(mean), 
        cv = as.vector(cv))
    for (i in c("q", "mean", "cv")) assign(i, arg.mat[, i])
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        p <- rep(NA, length(q))
    else {
        if (any(c(mean[!na.index], cv[!na.index]) < .Machine$double.eps)) 
            stop("All values of 'mean' and 'cv' must be positive.")
        sdlog <- sqrt(log(1 + cv^2))
        meanlog <- log(mean) - (sdlog^2)/2
        p <- plnorm(q = q, meanlog = meanlog, sdlog = sdlog, 
            lower.tail = lower.tail, log.p = log.p)
    }
    if (!is.null(names.q)) 
        names(p) <- rep(names.q, length = length(p))
    else names(p) <- NULL
    p
}
