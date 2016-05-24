predIntNormSimultaneousTestPower <-
function (n, df = n - 1, n.mean = 1, k = 1, m = 2, r = 1, rule = "k.of.m", 
    delta.over.sigma = 0, pi.type = "upper", conf.level = 0.95, 
    r.shifted = r, K.tol = .Machine$double.eps^0.5, integrate.args.list = NULL) 
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
            "and no Missing (NA), Infinite (Inf, -Inf),", "or Undefined (Nan) values."))
    if (!is.vector(n.mean, mode = "numeric") || !all(is.finite(n.mean)) || 
        any(n.mean < 1)) 
        stop(paste("'n.mean' must be a numeric vector", "with all elements greater than or equal to 1", 
            "and no Missing (NA), Infinite (Inf, -Inf),", "or Undefined (Nan) values."))
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
    if (!is.vector(delta.over.sigma, mode = "numeric") || any(is.na(delta.over.sigma))) 
        stop(paste("'delta.over.sigma' must be a numeric vector", 
            "with no Missing (NA) or Undefined (Nan) values."))
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
        n.mean = as.vector(n.mean), k = as.vector(k), m = as.vector(m), 
        r = as.vector(r), delta.over.sigma = as.vector(delta.over.sigma), 
        conf.level = as.vector(conf.level), r.shifted = as.vector(r.shifted))
    nrow.arg.mat <- nrow(arg.mat)
    length.rule <- length(rule)
    if (length.rule > nrow.arg.mat) 
        arg.mat <- arg.mat[rep(1:nrow.arg.mat, length.out = length.rule), 
            ]
    else rule <- rep(rule, length.out = nrow.arg.mat)
    for (i in c("n", "df", "n.mean", "k", "m", "r", "delta.over.sigma", 
        "conf.level", "r.shifted")) assign(i, arg.mat[, i])
    index <- rule == "k.of.m"
    if (any(index)) {
        if (any(k[index] > m[index])) 
            stop(paste("For cases where rule='k.of.m',", "all elements of 'k' must be less than or equal to", 
                "the corresponding elements of 'm'"))
    }
    index <- rule == "Modified.CA"
    m[index] <- 4
    N <- length(n)
    power <- numeric(N)
    index.0 <- delta.over.sigma == 0
    if (any(index.0)) 
        power[index.0] <- 1 - conf.level[index.0]
    if (!all(index.0)) {
        for (i in c("n", "df", "n.mean", "k", "m", "r", "delta.over.sigma", 
            "conf.level", "r.shifted")) assign(i, arg.mat[!index.0, 
            i])
        rule <- rule[!index.0]
        N <- length(n)
        power.sub <- numeric(N)
        index.easy <- rule == "k.of.m" & k == m & r == 1
        if (any(index.easy)) 
            power.sub[index.easy] <- predIntNormTestPower(n = n[index.easy], 
                df = df[index.easy], n.mean = n.mean[index.easy], 
                k = k[index.easy], delta.over.sigma = delta.over.sigma[index.easy], 
                pi.type = pi.type, conf.level = conf.level[index.easy])
        if (!all(index.easy)) {
            for (i in c("n", "df", "n.mean", "k", "m", "r", "delta.over.sigma", 
                "conf.level", "r.shifted")) assign(i, arg.mat[!index.0, 
                , drop = FALSE][!index.easy, i])
            rule <- rule[!index.easy]
            N <- length(n)
            if (all(n == n[1]) & all(df == df[1]) & all(k == 
                k[1]) & all(m == m[1]) & all(n.mean == n.mean[1]) & 
                all(r == r[1]) & all(rule == rule[1]) & all(conf.level == 
                conf.level[1])) {
                K <- predIntNormSimultaneousK(n = n[1], df = df[1], 
                  n.mean = n.mean[1], k = k[1], m = m[1], r = r[1], 
                  rule = rule[1], delta.over.sigma = 0, pi.type = pi.type, 
                  conf.level = conf.level[1], K.tol = K.tol)
                K <- rep(K, N)
            }
            else {
                K <- predIntNormSimultaneousK(n = n, df = df, 
                  n.mean = n.mean, k = k, m = m, r = r, rule = rule, 
                  delta.over.sigma = 0, pi.type = pi.type, conf.level = conf.level, 
                  K.tol = K.tol)
            }
            for (i in 1:N) power.sub[!index.easy][i] <- predIntNormSimultaneousTestPowerScalar(n = n[i], 
                df = df[i], n.mean = n.mean[i], K = K[i], k = k[i], 
                m = m[i], r = r.shifted[i], rule = rule[i], delta.over.sigma = delta.over.sigma[i], 
                pi.type = pi.type, conf.level = conf.level[i], 
                integrate.args.list = integrate.args.list)
        }
        power[!index.0] <- power.sub
    }
    power
}
