predIntNormHalfWidth <-
function (n, df = n - 1, n.mean = 1, k = 1, sigma.hat = 1, method = "Bonferroni", 
    conf.level = 0.95) 
{
    method <- match.arg(method, c("Bonferroni", "exact"))
    if (!is.vector(n, mode = "numeric") || any(is.na(n)) || any(n < 
        2)) 
        stop(paste("'n' must be a numeric vector", "with all elements greater than or equal to 2", 
            "and no Missing (NA), or Undefined (Nan) values."))
    if (!is.vector(df, mode = "numeric") || any(is.na(df)) || 
        any(df < 1)) 
        stop(paste("'df' must be a numeric vector", "with all elements greater than or equal to 1", 
            "and no Missing (NA), or Undefined (Nan) values."))
    if (!is.vector(n.mean, mode = "numeric") || !all(is.finite(n.mean)) || 
        !all(n.mean == trunc(n.mean)) || any(n.mean < 1)) 
        stop(paste("'n.mean' must be a numeric vector of integers", 
            "with all elements greater than 0", "and no Missing (NA), Infinite(-Inf, Inf),", 
            "or Undefined (Nan) values."))
    if (!is.vector(k, mode = "numeric") || !all(is.finite(k)) || 
        !all(k == trunc(k)) || any(k < 1)) 
        stop(paste("'k' must be a numeric vector of integers", 
            "with all elements greater than 0", "and no Missing (NA), Infinite(-Inf, Inf),", 
            "or Undefined (Nan) values."))
    if (!is.vector(sigma.hat, mode = "numeric") || !all(is.finite(sigma.hat)) || 
        any(sigma.hat < .Machine$double.eps)) 
        stop(paste("'sigma.hat' must be a numeric vector", "with all positive elements", 
            "and no Missing (NA), Infinite(-Inf, Inf),", "or Undefined (Nan) values."))
    if (!is.vector(conf.level, mode = "numeric") || !all(is.finite(conf.level)) || 
        any(conf.level <= .Machine$double.eps) || any(conf.level >= 
        1 - .Machine$double.eps)) 
        stop(paste("'conf.level' must be a numeric vector", "with all elements between 0 and 1", 
            "and no Missing (NA), Infinite(-Inf, Inf),", "or Undefined (Nan) values."))
    arg.mat <- cbind.no.warn(n = as.vector(n), df = as.vector(df), 
        n.mean = as.vector(n.mean), k = as.vector(k), sigma.hat = as.vector(sigma.hat), 
        conf.level = as.vector(conf.level))
    for (i in c("n", "df", "n.mean", "k", "sigma.hat", "conf.level")) assign(i, 
        arg.mat[, i])
    N <- length(n)
    K <- numeric(N)
    for (i in 1:N) K[i] <- predIntNormK(n = n[i], df = df[i], 
        n.mean = n.mean[i], k = k[i], method = method, pi.type = "two-sided", 
        conf.level = conf.level[i])
    half.width <- K * sigma.hat
    names(half.width) <- NULL
    half.width
}
