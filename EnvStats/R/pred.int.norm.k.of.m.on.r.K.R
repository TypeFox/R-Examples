pred.int.norm.k.of.m.on.r.K <-
function (n, df = n - 1, n.mean = 1, k = 1, m = 1, r = 1, delta.over.sigma = 0, 
    pi.type = c("upper", "lower"), conf.level = 0.95, K.tol = .Machine$double.eps^(1/2), 
    integrate.args.list = NULL) 
{
    if (!is.vector(n, mode = "numeric") || length(n) != 1 || 
        !is.vector(df, mode = "numeric") || length(df) != 1 || 
        !is.vector(k, mode = "numeric") || length(k) != 1 || 
        !is.vector(m, mode = "numeric") || length(m) != 1 || 
        !is.vector(n.mean, mode = "numeric") || length(n.mean) != 
        1 || !is.vector(r, mode = "numeric") || length(r) != 
        1 || !is.vector(delta.over.sigma, mode = "numeric") || 
        length(delta.over.sigma) != 1 || !is.vector(conf.level, 
        mode = "numeric") || length(conf.level) != 1) 
        stop(paste("'n', 'df', 'k', 'm', 'n.mean', 'r', 'delta.over.sigma',", 
            "and 'conf.level' must be numeric scalars"))
    if (n < 2 || n != trunc(n)) 
        stop("'n' must be an integer greater than or equal to 2")
    if (df < 1 || df != trunc(df)) 
        stop("'df' must be an integer greater than or equal to 1")
    if (k < 1) 
        stop("'k' must be greater than or equal to 1")
    if (m < 1) 
        stop("'m' must be greater than or equal to 1")
    if (k > m) 
        stop("'k' must be between 1 and 'm'")
    if (n.mean < 1) 
        stop("'n.mean' must be greater than or equal to 1")
    if (r < 1) 
        stop("'r' must be greater than or equal to 1")
    if (!is.finite(delta.over.sigma)) 
        stop("'delta.over.sigma' must be finite")
    if (conf.level <= 0 || conf.level >= 1) 
        stop("'conf.level' must be between 0 and 1")
    pi.type <- match.arg(pi.type)
    if (k == m & r == 1 & delta.over.sigma == 0) {
        K <- predIntNormK(n = n, df = df, n.mean = n.mean, k = m, 
            method = "exact", pi.type = pi.type, conf.level = conf.level)
    }
    else {
        fcn.to.min <- function(K, n.weird, df.weird, n.mean, 
            k.weird, m.weird, r.weird, delta.over.sigma, conf.level, 
            integrate.args.list) {
            (conf.level - pred.int.norm.k.of.m.on.r.prob(n = n.weird, 
                df = df.weird, n.mean = n.mean, K = K, delta.over.sigma = delta.over.sigma, 
                k = k.weird, m = m.weird, r = r.weird, integrate.args.list = integrate.args.list))^2
        }
        K <- nlminb(start = 1, objective = fcn.to.min, lower = 0, 
            control = list(x.tol = K.tol), n.weird = n, df.weird = df, 
            n.mean = n.mean, k.weird = k, m.weird = m, r.weird = r, 
            delta.over.sigma = delta.over.sigma, conf.level = conf.level, 
            integrate.args.list = integrate.args.list)$par
    }
    K
}
