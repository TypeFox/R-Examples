predIntLnormAltTestPower <-
function (n, df = n - 1, n.geomean = 1, k = 1, ratio.of.means = 1, 
    cv = 1, pi.type = "upper", conf.level = 0.95) 
{
    pi.type <- match.arg(pi.type, c("upper", "lower"))
    if (!is.vector(n, mode = "numeric") || !all(is.finite(n)) || 
        any(n < 2)) 
        stop(paste("'n' must be a numeric vector", "with all elements greater than or equal to 2", 
            "and no Missing (NA), Infinite (-Inf, Inf),", "or Undefined (Nan) values."))
    if (!is.vector(df, mode = "numeric") || !all(is.finite(df)) || 
        any(df < 1)) 
        stop(paste("'df' must be a numeric vector", "with all elements greater than or equal to 1", 
            "and no Missing (NA), Infinite (-Inf, Inf),", "or Undefined (Nan) values."))
    if (!is.vector(n.geomean, mode = "numeric") || !all(is.finite(n.geomean)) || 
        any(n.geomean < 1)) 
        stop(paste("'n.geomean' must be a numeric vector", "with all elements greater than or equal to 1", 
            "and no Missing (NA), Infinite (-Inf, Inf),", "or Undefined (Nan) values."))
    if (!is.vector(k, mode = "numeric") || !all(is.finite(k)) || 
        !all(k == trunc(k)) || any(k < 1)) 
        stop(paste("'k' must be a numeric vector of integers", 
            "with all elements greater than 0", "and no Missing (NA), Infinite(-Inf, Inf),", 
            "or Undefined (Nan) values."))
    if (!is.vector(ratio.of.means, mode = "numeric") || !all(is.finite(ratio.of.means)) || 
        any(ratio.of.means <= .Machine$double.eps)) 
        stop(paste("'ratio.of.means' must be a numeric vector", 
            "of positive values", "with no Missing (NA), Infinite(-Inf, Inf),", 
            "or Undefined (Nan) values."))
    if (!is.vector(cv, mode = "numeric") || !all(is.finite(cv)) || 
        any(cv <= .Machine$double.eps)) 
        stop(paste("'cv' must be a numeric vector", "of positive values", 
            "with no Missing (NA), Infinite(-Inf, Inf),", "or Undefined (Nan) values."))
    if (!is.vector(conf.level, mode = "numeric") || !all(is.finite(conf.level)) || 
        any(conf.level <= .Machine$double.eps) || any(conf.level >= 
        1 - .Machine$double.eps)) 
        stop(paste("'conf.level' must be a numeric vector", "with all elements between 0 and 1", 
            "and no Missing (NA), Infinite(-Inf, Inf),", "or Undefined (Nan) values."))
    arg.mat <- cbind.no.warn(n = as.vector(n), df = as.vector(df), 
        n.geomean = as.vector(n.geomean), k = as.vector(k), ratio.of.means = as.vector(ratio.of.means), 
        cv = as.vector(cv), conf.level = as.vector(conf.level))
    for (i in c("n", "df", "n.geomean", "k", "ratio.of.means", 
        "cv", "conf.level")) assign(i, arg.mat[, i])
    delta.over.sigma <- log(ratio.of.means)/sqrt(log(cv^2 + 1))
    predIntNormTestPower(n = n, df = df, n.mean = n.geomean, 
        k = k, delta.over.sigma = delta.over.sigma, pi.type = pi.type, 
        conf.level = conf.level)
}
