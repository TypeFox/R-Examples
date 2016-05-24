predIntNormSimultaneousTestPowerScalar <-
function (n, df = n - 1, n.mean = 1, K, k = 1, m = 2, r = 1, 
    rule = c("k.of.m", "CA", "Modified.CA"), delta.over.sigma = 0, 
    pi.type = c("upper", "lower"), conf.level = 0.95, integrate.args.list = NULL) 
{
    if (!is.vector(K, mode = "numeric") || length(K) != 1 || 
        K <= 0) 
        stop("'K' must be a positive scalar")
    rule <- match.arg(rule)
    switch(rule, k.of.m = {
        if (!is.vector(n, mode = "numeric") || length(n) != 1 || 
            !is.vector(df, mode = "numeric") || length(df) != 
            1 || !is.vector(k, mode = "numeric") || length(k) != 
            1 || !is.vector(m, mode = "numeric") || length(m) != 
            1 || !is.vector(n.mean, mode = "numeric") || length(n.mean) != 
            1 || !is.vector(r, mode = "numeric") || length(r) != 
            1 || !is.vector(delta.over.sigma, mode = "numeric") || 
            length(delta.over.sigma) != 1 || !is.vector(conf.level, 
            mode = "numeric") || length(conf.level) != 1) stop(paste("'n', 'df', 'k', 'm', 'n.mean', 'r', 'delta.over.sigma',", 
            "and 'conf.level' must be numeric scalars"))
    }, CA = {
        if (!is.vector(n, mode = "numeric") || length(n) != 1 || 
            !is.vector(df, mode = "numeric") || length(df) != 
            1 || !is.vector(m, mode = "numeric") || length(m) != 
            1 || !is.vector(n.mean, mode = "numeric") || length(n.mean) != 
            1 || !is.vector(r, mode = "numeric") || length(r) != 
            1 || !is.vector(delta.over.sigma, mode = "numeric") || 
            length(delta.over.sigma) != 1 || !is.vector(conf.level, 
            mode = "numeric") || length(conf.level) != 1) stop(paste("'n', 'df', 'm', 'n.mean', 'r', 'delta.over.sigma',", 
            "and 'conf.level' must be numeric scalars"))
    }, Modified.CA = {
        if (!is.vector(n, mode = "numeric") || length(n) != 1 || 
            !is.vector(df, mode = "numeric") || length(df) != 
            1 || !is.vector(m, mode = "numeric") || length(m) != 
            1 || !is.vector(n.mean, mode = "numeric") || length(n.mean) != 
            1 || !is.vector(r, mode = "numeric") || length(r) != 
            1 || !is.vector(delta.over.sigma, mode = "numeric") || 
            length(delta.over.sigma) != 1 || !is.vector(conf.level, 
            mode = "numeric") || length(conf.level) != 1) stop(paste("'n', 'df', 'm', 'n.mean', 'r', 'delta.over.sigma',", 
            "and 'conf.level' must be numeric scalars"))
        m <- 4
    })
    if (!is.finite(n) || n < 2) 
        stop("'n' must be finite and greater than or equal to 2")
    if (!is.finite(df) || df < 1) 
        stop("'df' must be finite and greater than or equal to 1")
    if (!is.finite(m) || m < 1) 
        stop("'m' must be finite and greater than or equal to 1")
    if (!is.finite(n.mean) || n.mean < 1) 
        stop("'n.mean' must be finite and greater than or equal to 1")
    if (rule == "k.of.m") {
        if (k < 1) 
            stop("'k' must be finite and greater than or equal to 1")
        if (k > m) 
            stop("'k' must be between 1 and 'm'")
    }
    if (!is.finite(r) || r < 1) 
        stop("'r' must be finite and greater than or equal to 1")
    if (is.na(delta.over.sigma)) 
        stop(paste("'delta.over.sigma' cannot be a missing (NA)", 
            "or undefined (NaN) value"))
    if (conf.level <= 0 || conf.level >= 1) 
        stop("'conf.level' must be between 0 and 1")
    pi.type <- match.arg(pi.type)
    if (delta.over.sigma == 0) 
        power <- 1 - conf.level
    else if (!is.finite(delta.over.sigma)) {
        if (pi.type == "upper") {
            if (delta.over.sigma == Inf) 
                power <- 1
            else power <- 0
        }
        else {
            if (delta.over.sigma == -Inf) 
                power <- 1
            else power <- 0
        }
    }
    else {
        if (rule == "k.of.m" & k == m & r == 1) {
            power <- predIntNormTestPowerScalar(n = n, df = df, 
                n.mean = n.mean, K = K, k = k, delta.over.sigma = delta.over.sigma, 
                pi.type = pi.type, conf.level = conf.level)
        }
        else {
            power <- switch(pi.type, upper = 1 - predIntNormSimultaneousProb(n = n, 
                df = df, n.mean = n.mean, K = K, delta.over.sigma = delta.over.sigma, 
                k = k, m = m, r = r, rule = rule, integrate.args.list = integrate.args.list), 
                lower = 1 - predIntNormSimultaneousProb(n = n, 
                  df = df, n.mean = n.mean, K = K, delta.over.sigma = -delta.over.sigma, 
                  k = k, m = m, r = r, rule = rule, integrate.args.list = integrate.args.list))
        }
    }
    power
}
