predIntNparSimultaneousN <-
function (n.median = 1, k = 1, m = 2, r = 1, rule = "k.of.m", 
    lpl.rank = ifelse(pi.type == "upper", 0, 1), n.plus.one.minus.upl.rank = ifelse(pi.type == 
        "lower", 0, 1), pi.type = "upper", conf.level = 0.95, 
    n.max = 5000, integrate.args.list = NULL, maxiter = 1000) 
{
    rule <- match.arg(rule, c("k.of.m", "CA", "Modified.CA"), 
        several.ok = TRUE)
    if (!is.vector(n.median, mode = "numeric") || length(n.median) != 
        1 || !is.finite(n.median) || n.median != trunc(n.median) || 
        n.median < 1 || !is.odd(n.median)) 
        stop("'n.median' must be a positive odd integer")
    if (!is.vector(n.median, mode = "numeric") || !all(is.finite(n.median)) || 
        any(n.median < 1) || !all(n.median == trunc(n.median)) || 
        !all(is.odd(n.median))) 
        stop("'n.median' must be a numeric vector of positive odd integers")
    if (!is.vector(k, mode = "numeric") || !all(is.finite(k)) || 
        any(k < 1)) 
        stop(paste("'k' must be a numeric vector", "with all elements greater than or equal to 1", 
            "and no Missing (NA), Infinite (-Inf, Inf),", "or Undefined (Nan) values."))
    if (!is.vector(m, mode = "numeric") || !all(is.finite(m)) || 
        any(m < 1)) 
        stop(paste("'m' must be a numeric vector", "with all elements greater than or equal to 1", 
            "and no Missing (NA), Infinite (Inf, -Inf),", "or Undefined (Nan) values."))
    if (!is.vector(r, mode = "numeric") || !all(is.finite(r)) || 
        any(r < 1)) 
        stop(paste("'r' must be a numeric vector", "with all elements greater than or equal to 1", 
            "and no Missing (NA), Infinite (-Inf, Inf),", "or Undefined (Nan) values."))
    if (!is.vector(r, mode = "numeric") || !all(is.finite(r)) || 
        any(r < 1)) 
        stop(paste("'r' must be a numeric vector", "with all elements greater than or equal to 1", 
            "and no Missing (NA), Infinite (-Inf, Inf),", "or Undefined (Nan) values."))
    if (!is.vector(conf.level, mode = "numeric") || !all(is.finite(conf.level)) || 
        any(conf.level <= .Machine$double.eps) || any(conf.level >= 
        1 - .Machine$double.eps)) 
        stop(paste("'conf.level' must be a numeric vector", "with all elements between 0 and 1", 
            "and no Missing (NA), Infinite(-Inf, Inf),", "or Undefined (Nan) values."))
    if (!is.vector(n.max, mode = "numeric") || length(n.max) != 
        1 || !is.finite(n.max) || n.max != trunc(n.max) || n.max < 
        2) 
        stop("'n.max' must be a positive integer greater than 1")
    if (!is.vector(maxiter, mode = "numeric") || length(maxiter) != 
        1 || !is.finite(maxiter) || maxiter != trunc(maxiter) || 
        maxiter < 2) 
        stop("'maxiter' must be a positive integer greater than 1")
    pi.type <- match.arg(pi.type, c("upper", "lower"))
    if (pi.type == "upper") 
        lpl.rank <- 0
    else n.plus.one.minus.upl.rank <- 0
    if (!is.vector(lpl.rank, mode = "numeric") || !all(is.finite(lpl.rank)) || 
        any(lpl.rank != trunc(lpl.rank)) || any(lpl.rank < 0 | 
        lpl.rank >= n.max)) 
        stop(paste("'lpl.rank' must be a vector of non-negative integers", 
            "with all elements less than 'n.max'"))
    if (pi.type == "lower" & any(lpl.rank < 1)) 
        stop("When pi.type='lower', all values of 'lpl.rank' must be positive integers")
    if (!is.vector(n.plus.one.minus.upl.rank, mode = "numeric") || 
        !all(is.finite(n.plus.one.minus.upl.rank)) || any(n.plus.one.minus.upl.rank != 
        trunc(n.plus.one.minus.upl.rank)) || any(n.plus.one.minus.upl.rank < 
        0 | n.plus.one.minus.upl.rank >= n.max)) 
        stop(paste("'n.plus.one.minus.upl.rank' must be a vector of non-negative", 
            "integers less than 'n.max'"))
    if (pi.type == "upper" & any(n.plus.one.minus.upl.rank < 
        1)) 
        stop("When pi.type='upper' all values of 'n.plus.one.minus.upl.rank' must be positive integers")
    arg.mat <- cbind.no.warn(n.median = as.vector(n.median), 
        k = as.vector(k), m = as.vector(m), r = as.vector(r), 
        lpl.rank = as.vector(lpl.rank), n.plus.one.minus.upl.rank = as.vector(n.plus.one.minus.upl.rank), 
        conf.level = as.vector(conf.level))
    nrow.arg.mat <- nrow(arg.mat)
    length.rule <- length(rule)
    if (length.rule > nrow.arg.mat) 
        arg.mat <- arg.mat[rep(1:nrow.arg.mat, length.out = length.rule), 
            ]
    else rule <- rep(rule, length.out = nrow.arg.mat)
    for (i in c("n.median", "k", "m", "r", "lpl.rank", "n.plus.one.minus.upl.rank", 
        "conf.level")) assign(i, arg.mat[, i])
    index <- rule == "k.of.m"
    if (any(index)) {
        if (any(k[index] > m[index])) 
            stop(paste("For cases where rule='k.of.m',", "all elements of 'k' must be less than or equal to", 
                "the corresponding elements of 'm'"))
        if (any(k[index] > 1)) 
            stop(paste("Determining sample size for the k-of-m rule", 
                "with k > 1 is not yet implemented,", "so currently all values of 'k' must be equal to 1"))
    }
    index <- rule == "Modified.CA"
    m[index] <- 4
    N <- length(k)
    n <- numeric(N)
    for (i in 1:N) {
        n[i] <- predIntNparSimultaneousNScalar(n.median = n.median[i], 
            k = k[i], m = m[i], r = r[i], rule = rule[i], lpl.rank = lpl.rank[i], 
            n.plus.one.minus.upl.rank = n.plus.one.minus.upl.rank[i], 
            pi.type = pi.type, conf.level = conf.level[i], n.max = n.max, 
            integrate.args.list = integrate.args.list, maxiter = maxiter)
    }
    n
}
