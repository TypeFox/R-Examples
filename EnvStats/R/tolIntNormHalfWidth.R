tolIntNormHalfWidth <-
function (n, sigma.hat = 1, coverage = 0.95, cov.type = "content", 
    conf.level = 0.95, method = "wald.wolfowitz") 
{
    cov.type <- match.arg(cov.type, c("content", "expectation"))
    method <- match.arg(method, c("exact", "wald.wolfowitz"))
    if (cov.type == "content") {
        if (!is.vector(n, mode = "numeric") || !all(is.finite(n)) || 
            any(n < 2)) 
            stop(paste("'n' must be a numeric vector", "with all elements greater than 1", 
                "and no Missing (NA), Infinite(-Inf, Inf),", 
                "or Undefined (Nan) values."))
    }
    else {
        if (!is.vector(n, mode = "numeric") || any(is.na(n)) || 
            any(n < 2)) 
            stop(paste("'n' must be a numeric vector", "with all elements greater than 1", 
                "and no Missing (NA), or Undefined (Nan) values."))
    }
    if (!is.vector(sigma.hat, mode = "numeric") || !all(is.finite(sigma.hat)) || 
        any(sigma.hat < .Machine$double.eps)) 
        stop(paste("'sigma.hat' must be a numeric vector", "with all positive elements", 
            "and no Missing (NA), Infinite(-Inf, Inf),", "or Undefined (Nan) values."))
    if (!is.vector(coverage, mode = "numeric") || !all(is.finite(coverage)) || 
        any(coverage <= .Machine$double.eps) || any(coverage >= 
        1 - .Machine$double.eps)) 
        stop(paste("'coverage' must be a numeric vector", "with all elements between 0 and 1", 
            "and no Missing (NA), Infinite(-Inf, Inf),", "or Undefined (Nan) values."))
    if (cov.type == "content") {
        if (!is.vector(conf.level, mode = "numeric") || !all(is.finite(conf.level)) || 
            any(conf.level <= .Machine$double.eps) || any(conf.level >= 
            1 - .Machine$double.eps)) 
            stop(paste("'conf.level' must be a numeric vector", 
                "with all elements between 0 and 1", "and no Missing (NA), Infinite(-Inf, Inf),", 
                "or Undefined (Nan) values."))
        arg.mat <- cbind.no.warn(n = as.vector(n), sigma.hat = as.vector(sigma.hat), 
            coverage = as.vector(coverage), conf.level = as.vector(conf.level))
        for (i in c("n", "sigma.hat", "coverage", "conf.level")) assign(i, 
            arg.mat[, i])
        N <- length(n)
        K <- numeric(N)
        for (i in 1:N) K[i] <- tolIntNormK(n = n[i], coverage = coverage[i], 
            cov.type = cov.type, ti.type = "two-sided", conf.level = conf.level[i], 
            method = method)
    }
    else {
        arg.mat <- cbind.no.warn(n = as.vector(n), sigma.hat = as.vector(sigma.hat), 
            coverage = as.vector(coverage))
        for (i in c("n", "sigma.hat", "coverage")) assign(i, 
            arg.mat[, i])
        N <- length(n)
        K <- numeric(N)
        for (i in 1:N) K[i] <- tolIntNormK(n = n[i], coverage = coverage[i], 
            cov.type = cov.type, ti.type = "two-sided", method = method)
    }
    half.width <- K * sigma.hat
    names(half.width) <- NULL
    half.width
}
