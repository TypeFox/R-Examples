predIntNparN <-
function (k = m, m = 1, lpl.rank = ifelse(pi.type == "upper", 
    0, 1), n.plus.one.minus.upl.rank = ifelse(pi.type == "lower", 
    0, 1), pi.type = "two.sided", conf.level = 0.95, n.max = 5000, 
    maxiter = 1000) 
{
    if (!is.vector(m, mode = "numeric") || !all(is.finite(m)) || 
        any(m != trunc(m)) || any(m < 1)) 
        stop("'m' must be a vector of positive integers")
    if (!is.vector(k, mode = "numeric") || !all(is.finite(k)) || 
        any(k != trunc(k)) || any(k < 1) || any(k > m)) 
        stop(paste("'k' must be a vector of positive integers,", 
            "and all values of 'k' must be between", "1 and the corresponding value of 'm'"))
    pi.type <- match.arg(pi.type, c("two.sided", "lower", "upper"))
    if (pi.type == "upper") 
        lpl.rank <- 0
    else if (pi.type == "lower") 
        n.plus.one.minus.upl.rank <- 0
    if (!is.vector(lpl.rank, mode = "numeric") || !all(is.finite(lpl.rank)) || 
        any(lpl.rank != trunc(lpl.rank)) || any(lpl.rank < 0)) 
        stop("'lpl.rank' must be a vector of non-negative integers")
    if (pi.type %in% c("two.sided", "lower") & any(lpl.rank < 
        1)) 
        stop("When pi.type='two.sided' or pi.type='lower', all values of 'lpl.rank' must be positive integers")
    if (!is.vector(n.plus.one.minus.upl.rank, mode = "numeric") || 
        !all(is.finite(n.plus.one.minus.upl.rank)) || any(n.plus.one.minus.upl.rank != 
        trunc(n.plus.one.minus.upl.rank)) || any(n.plus.one.minus.upl.rank < 
        0)) 
        stop("'n.plus.one.minus.upl.rank' must be a vector of non-negative integers")
    if (pi.type %in% c("two.sided", "upper") & any(n.plus.one.minus.upl.rank < 
        1)) 
        stop("When pi.type='two.sided' or pi.type='upper' all values of 'n.plus.one.minus.upl.rank' must be positive integers")
    if (any(conf.level <= .Machine$double.eps) || any(conf.level >= 
        1 - .Machine$double.eps)) 
        stop("All values of 'conf.level' must be greater than 0 and less than 1")
    if (!is.vector(maxiter, mode = "numeric") || length(maxiter) != 
        1 || !is.finite(maxiter) || maxiter != trunc(maxiter) || 
        maxiter < 1) 
        stop("'maxiter' must be a positive integer")
    if (!is.vector(n.max, mode = "numeric") || length(n.max) != 
        1 || !is.finite(n.max) || n.max != trunc(n.max) || n.max < 
        2) 
        stop("'n.max' must be a positive integer greater than 1")
    if (!is.vector(maxiter, mode = "numeric") || length(maxiter) != 
        1 || !is.finite(maxiter) || maxiter != trunc(maxiter) || 
        maxiter < 2) 
        stop("'maxiter' must be a positive integer greater than 1")
    arg.mat <- cbind.no.warn(k = as.vector(k), m = as.vector(m), 
        lpl.rank = as.vector(lpl.rank), n.plus.one.minus.upl.rank = as.vector(n.plus.one.minus.upl.rank), 
        conf.level = as.vector(conf.level))
    for (i in c("k", "m", "lpl.rank", "n.plus.one.minus.upl.rank", 
        "conf.level")) assign(i, arg.mat[, i])
    N <- length(k)
    n.vec <- numeric(N)
    for (i in 1:N) {
        n.vec[i] <- predIntNparNScalar(k = k[i], m = m[i], lpl.rank = lpl.rank[i], 
            n.plus.one.minus.upl.rank = n.plus.one.minus.upl.rank[i], 
            pi.type = pi.type, conf.level = conf.level[i], n.max = n.max, 
            maxiter = maxiter)
    }
    names(n.vec) <- NULL
    n.vec
}
