qlnormAlt <-
function (p, mean = exp(1/2), cv = sqrt(exp(1) - 1), lower.tail = TRUE, 
    log.p = FALSE) 
{
    names.p <- names(p)
    arg.mat <- cbind.no.warn(p = as.vector(p), mean = as.vector(mean), 
        cv = as.vector(cv))
    for (i in c("p", "mean", "cv")) assign(i, arg.mat[, i])
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        q <- rep(NA, length(p))
    else {
        if (any(c(mean[!na.index], cv[!na.index]) < .Machine$double.eps)) 
            stop("All values of 'mean' and 'cv' must be positive.")
        sdlog <- sqrt(log(1 + cv^2))
        meanlog <- log(mean) - (sdlog^2)/2
        q <- qlnorm(p = p, meanlog = meanlog, sdlog = sdlog, 
            lower.tail = lower.tail, log.p = log.p)
    }
    if (!is.null(names.p)) 
        names(q) <- rep(names.p, length = length(q))
    else names(q) <- NULL
    q
}
