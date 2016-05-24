ciNparN <-
function (p = 0.5, lcl.rank = ifelse(ci.type == "upper", 0, 1), 
    n.plus.one.minus.ucl.rank = ifelse(ci.type == "lower", 0, 
        1), ci.type = "two.sided", conf.level = 0.95) 
{
    if (!is.vector(p, mode = "numeric") || !all(is.finite(p)) || 
        any(p <= 0) || any(p >= 1)) 
        stop("All values of 'p' must be greater than 0 and less than 1.")
    ci.type <- match.arg(ci.type, c("two.sided", "lower", "upper"))
    if (!is.vector(lcl.rank, mode = "numeric") || !all(is.finite(lcl.rank)) || 
        any(lcl.rank != trunc(lcl.rank)) || any(lcl.rank < 0)) 
        stop("All values of 'lcl.rank' must be non-negative integers")
    if (!is.vector(n.plus.one.minus.ucl.rank, mode = "numeric") || 
        !all(is.finite(n.plus.one.minus.ucl.rank)) || any(n.plus.one.minus.ucl.rank != 
        trunc(n.plus.one.minus.ucl.rank)) || any(n.plus.one.minus.ucl.rank < 
        0)) 
        stop("All values of 'n.plus.one.minus.ucl.rank' must be non-negative integers")
    if (!is.vector(conf.level, mode = "numeric") || !all(is.finite(conf.level)) || 
        any(conf.level <= 0) || any(conf.level >= 1)) 
        stop("All values of 'conf.level' must be greater than 0 and less than 1.")
    arg.mat <- cbind.no.warn(p = as.vector(p), lcl.rank = as.vector(lcl.rank), 
        n.plus.one.minus.ucl.rank = as.vector(n.plus.one.minus.ucl.rank), 
        conf.level = as.vector(conf.level))
    for (i in c("p", "lcl.rank", "n.plus.one.minus.ucl.rank", 
        "conf.level")) assign(i, arg.mat[, i])
    N <- length(p)
    n.vec <- numeric(N)
    for (i in 1:N) {
        n.vec[i] <- ciNparNScalar(p = p[i], lcl.rank = lcl.rank[i], 
            n.plus.one.minus.ucl.rank = n.plus.one.minus.ucl.rank[i], 
            ci.type = ci.type, conf.level = conf.level[i])
    }
    names(n.vec) <- NULL
    n.vec
}
