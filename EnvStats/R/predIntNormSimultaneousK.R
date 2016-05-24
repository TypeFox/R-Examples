predIntNormSimultaneousK <-
function (n, df = n - 1, n.mean = 1, k = 1, m = 2, r = 1, rule = "k.of.m", 
    delta.over.sigma = 0, pi.type = "upper", conf.level = 0.95, 
    K.tol = .Machine$double.eps^0.5, integrate.args.list = NULL) 
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
    if (!is.vector(k, mode = "numeric") || !all(is.finite(k)) || 
        any(k < 1)) 
        stop(paste("'k' must be a numeric vector", "with all elements greater than or eqal to 1", 
            "and no Missing (NA), Infinite (-Inf, Inf),", "or Undefined (Nan) values."))
    if (!is.vector(m, mode = "numeric") || !all(is.finite(m)) || 
        any(m < 1)) 
        stop(paste("'m' must be a numeric vector", "with all elements greater than or equal to 1", 
            "and no Missing (NA), Infinite (Inf, -Inf),", "or Undefined (Nan) values."))
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
    arg.mat <- cbind.no.warn(n = as.vector(n), df = as.vector(df), 
        n.mean = as.vector(n.mean), k = as.vector(k), m = as.vector(m), 
        r = as.vector(r), delta.over.sigma = as.vector(delta.over.sigma), 
        conf.level = as.vector(conf.level))
    nrow.arg.mat <- nrow(arg.mat)
    length.rule <- length(rule)
    if (length.rule > nrow.arg.mat) 
        arg.mat <- arg.mat[rep(1:nrow.arg.mat, length.out = length.rule), 
            ]
    else rule <- rep(rule, length.out = nrow.arg.mat)
    for (i in c("n", "df", "n.mean", "k", "m", "r", "delta.over.sigma", 
        "conf.level")) assign(i, arg.mat[, i])
    index <- rule == "k.of.m"
    if (any(index)) {
        if (any(k[index] > m[index])) 
            stop(paste("For cases where rule='k.of.m',", "all elements of 'k' must be less than or equal to", 
                "the corresponding elements of 'm'"))
    }
    index <- rule == "Modified.CA"
    m[index] <- 4
    N <- length(n)
    K <- numeric(N)
    for (i in 1:N) {
        K[i] <- switch(rule[i], k.of.m = {
            pred.int.norm.k.of.m.on.r.K(n = n[i], df = df[i], 
                n.mean = n.mean[i], k = k[i], m = m[i], r = r[i], 
                delta.over.sigma = delta.over.sigma[i], pi.type = pi.type, 
                conf.level = conf.level[i], K.tol = K.tol, integrate.args.list = integrate.args.list)
        }, CA = {
            pred.int.norm.CA.on.r.K(n = n[i], df = df[i], n.mean = n.mean[i], 
                m = m[i], r = r[i], delta.over.sigma = delta.over.sigma[i], 
                pi.type = pi.type, conf.level = conf.level[i], 
                K.tol = K.tol)
        }, Modified.CA = {
            pred.int.norm.Modified.CA.on.r.K(n = n[i], df = df[i], 
                n.mean = n.mean[i], r = r[i], delta.over.sigma = delta.over.sigma[i], 
                pi.type = pi.type, conf.level = conf.level[i], 
                K.tol = K.tol)
        })
    }
    K
}
