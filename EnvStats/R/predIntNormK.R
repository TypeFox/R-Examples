predIntNormK <-
function (n, df = n - 1, n.mean = 1, k = 1, method = "Bonferroni", 
    pi.type = "two-sided", conf.level = 0.95) 
{
    if (!is.vector(n, mode = "numeric") || length(n) != 1 || 
        !is.vector(df, mode = "numeric") || length(df) != 1 || 
        !is.vector(k, mode = "numeric") || length(k) != 1 || 
        !is.vector(n.mean, mode = "numeric") || length(n.mean) != 
        1 || !is.vector(conf.level, mode = "numeric") || length(conf.level) != 
        1) 
        stop("'n', 'df', 'k', 'n.mean', and 'conf.level' must be numeric scalars")
    if (n < 2) 
        stop("'n' must be greater than or equal to 2")
    if (df < 1) 
        stop("'df' must be greater than or equal to 1")
    if (k < 1) 
        stop("'k' must be greater than or equal to 1")
    if (n.mean < 1 || n.mean != trunc(n.mean)) 
        stop("'n.mean' must be an integer greater than 0")
    if (conf.level <= 0 || conf.level >= 1) 
        stop("'conf.level' must be between 0 and 1")
    method <- match.arg(method, c("Bonferroni", "exact"))
    pi.type <- match.arg(pi.type, c("two-sided", "lower", "upper"))
    Dunnetts.K(n = n, df = df, k = k, m = n.mean, method = method, 
        type = pi.type, conf.level = conf.level)
}
