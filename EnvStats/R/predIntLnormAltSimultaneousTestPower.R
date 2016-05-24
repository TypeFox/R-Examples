predIntLnormAltSimultaneousTestPower <-
function (n, df = n - 1, n.geomean = 1, k = 1, m = 2, r = 1, 
    rule = "k.of.m", ratio.of.means = 1, cv = 1, pi.type = "upper", 
    conf.level = 0.95, r.shifted = r, K.tol = .Machine$double.eps^0.5, 
    integrate.args.list = NULL) 
{
    rule <- match.arg(rule, c("k.of.m", "CA", "Modified.CA"), 
        several.ok = TRUE)
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
    if (!is.vector(m, mode = "numeric") || !all(is.finite(m)) || 
        any(m < 1)) 
        stop(paste("'m' must be a numeric vector", "with all elements greater than or equal to 1", 
            "and no Missing (NA), Infinite (Inf, -Inf),", "or Undefined (Nan) values."))
    if (!is.vector(k, mode = "numeric") || !all(is.finite(k)) || 
        any(k < 1)) 
        stop(paste("'k' must be a numeric vector", "with all elements greater than or eqal to 1", 
            "and no Missing (NA), Infinite (-Inf, Inf),", "or Undefined (Nan) values."))
    if (!is.vector(r, mode = "numeric") || !all(is.finite(r)) || 
        any(r < 1)) 
        stop(paste("'r' must be a numeric vector", "with all elements greater than or equal to 1", 
            "and no Missing (NA), Infinite (-Inf, Inf),", "or Undefined (Nan) values."))
    if (!is.vector(ratio.of.means, mode = "numeric") || any(is.na(ratio.of.means)) || 
        any(ratio.of.means <= .Machine$double.eps)) 
        stop(paste("'ratio.of.means' must be a numeric vector", 
            "of positive values", "with no Missing (NA) or Undefined (Nan) values."))
    if (!is.vector(cv, mode = "numeric") || !all(is.finite(cv)) || 
        any(cv <= .Machine$double.eps)) 
        stop(paste("'cv' must be a numeric vector", "of positive values", 
            "with no Missing (NA), Infinite(-Inf, Inf),", "or Undefined (Nan) values."))
    if (!is.vector(conf.level, mode = "numeric") || !all(is.finite(conf.level)) || 
        any(conf.level <= .Machine$double.eps) || any(conf.level >= 
        1 - .Machine$double.eps)) 
        stop(paste("'conf.level' must be a numeric vector", "with all elements between 0 and 1", 
            "and no Missing (NA), Infinite(-Inf, Inf),", "or Undefined (Nan) values."))
    if (!is.vector(r.shifted, mode = "numeric") || !all(is.finite(r.shifted)) || 
        !all(r.shifted == trunc(r.shifted)) || any(r.shifted < 
        1) || any(r.shifted > r)) 
        stop(paste("'r.shifted' must be a numeric vector of positive integers", 
            "with all values must be less than or equal to", 
            "the corresponding values of 'r'"))
    arg.mat <- cbind.no.warn(n = as.vector(n), df = as.vector(df), 
        n.geomean = as.vector(n.geomean), k = as.vector(k), m = as.vector(m), 
        r = as.vector(r), ratio.of.means = as.vector(ratio.of.means), 
        cv = as.vector(cv), conf.level = as.vector(conf.level), 
        r.shifted = as.vector(r.shifted))
    nrow.arg.mat <- nrow(arg.mat)
    length.rule <- length(rule)
    if (length.rule > nrow.arg.mat) 
        arg.mat <- arg.mat[rep(1:nrow.arg.mat, length.out = length.rule), 
            ]
    else rule <- rep(rule, length.out = nrow.arg.mat)
    for (i in c("n", "df", "n.geomean", "k", "m", "r", "ratio.of.means", 
        "cv", "conf.level", "r.shifted")) assign(i, arg.mat[, 
        i])
    index <- rule == "k.of.m"
    if (any(index)) {
        if (any(k[index] > m[index])) 
            stop(paste("For cases where rule='k.of.m',", "all elements of 'k' must be less than or equal to", 
                "the corresponding elements of 'm'"))
    }
    index <- rule == "Modified.CA"
    m[index] <- 4
    delta.over.sigma <- log(ratio.of.means)/sqrt(log(cv^2 + 1))
    predIntNormSimultaneousTestPower(n = n, df = df, n.mean = n.geomean, 
        k = k, m = m, r = r, rule = rule, delta.over.sigma = delta.over.sigma, 
        pi.type = pi.type, conf.level = conf.level, r.shifted = r.shifted, 
        K.tol = K.tol, integrate.args.list = integrate.args.list)
}
