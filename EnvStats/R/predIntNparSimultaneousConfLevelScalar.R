predIntNparSimultaneousConfLevelScalar <-
function (n, n.median = 1, k = 1, m = 2, r = 1, rule = c("k.of.m", 
    "CA", "Modified.CA"), lpl.rank = ifelse(pi.type == "upper", 
    0, 1), n.plus.one.minus.upl.rank = ifelse(pi.type == "lower", 
    0, 1), pi.type = "upper", integrate.args.list = NULL) 
{
    rule <- match.arg(rule)
    pi.type <- match.arg(pi.type, c("upper", "lower"))
    if (!is.vector(n, mode = "numeric") || length(n) != 1 || 
        n < 2 || n != trunc(n)) 
        stop("'n' must be an integer greater than 1")
    if (!is.vector(n.median, mode = "numeric") || length(n.median) != 
        1 || n.median < 1 || n.median != trunc(n.median) || !is.odd(n.median)) 
        stop("'n.median' must be an odd positive integer")
    switch(rule, k.of.m = {
        if (!is.vector(k, mode = "numeric") || length(k) != 1 || 
            !is.vector(m, mode = "numeric") || length(m) != 1 || 
            !is.vector(r, mode = "numeric") || length(r) != 1) stop(paste("'k', 'm', and 'r'", 
            "must be numeric scalars"))
    }, CA = {
        if (!is.vector(m, mode = "numeric") || length(m) != 1 || 
            !is.vector(r, mode = "numeric") || length(r) != 1) stop(paste("'m', and 'r'", 
            "must be numeric scalars"))
    }, Modified.CA = {
        if (!is.vector(r, mode = "numeric") || length(r) != 1) stop("'r' must be a numeric scalar")
        m <- 4
    })
    if (m < 1 || m != trunc(m)) 
        stop("'m' must be a positive integer")
    if (rule == "k.of.m") {
        if (k < 1 || k != trunc(k)) 
            stop("'k' must be a positive integer")
        if (k > m) 
            stop("'k' must be between >=1 and <='m'")
    }
    if (r < 1 || r != trunc(r)) 
        stop("'r' must be a positive integer")
    pi.type <- match.arg(pi.type, c("upper", "lower"))
    if (pi.type == "upper") 
        lpl.rank <- 0
    else n.plus.one.minus.upl.rank <- 0
    if (!is.numeric(lpl.rank) || length(lpl.rank) > 1 || lpl.rank != 
        trunc(lpl.rank) || lpl.rank < 0 || lpl.rank >= n) 
        stop("'lpl.rank' must be a non-negative integer less than 'n'")
    if (pi.type == "lower" & lpl.rank < 1) 
        stop("When pi.type='lower', 'lpl.rank' must be a positive integer")
    if (!is.numeric(n.plus.one.minus.upl.rank) || length(n.plus.one.minus.upl.rank) > 
        1 || n.plus.one.minus.upl.rank != trunc(n.plus.one.minus.upl.rank) || 
        n.plus.one.minus.upl.rank < 0 || n.plus.one.minus.upl.rank >= 
        n) 
        stop("'n.plus.one.minus.upl.rank' must be a non-negative integer less than 'n'")
    if (pi.type == "upper" & n.plus.one.minus.upl.rank < 1) 
        stop("When pi.type='upper', 'n.plus.one.minus.upl.rank' must be a positive integer")
    conf.level <- switch(rule, k.of.m = {
        pred.int.npar.k.of.m.on.r.conf.level(n = n, n.median = n.median, 
            k = k, m = m, r = r, lpl.rank = lpl.rank, upl.rank = n + 
                1 - n.plus.one.minus.upl.rank, pi.type = pi.type, 
            integrate.args.list = integrate.args.list)
    }, CA = {
        pred.int.npar.CA.on.r.conf.level(n = n, n.median = n.median, 
            m = m, r = r, lpl.rank = lpl.rank, upl.rank = n + 
                1 - n.plus.one.minus.upl.rank, pi.type = pi.type, 
            integrate.args.list = integrate.args.list)
    }, Modified.CA = {
        pred.int.npar.Modified.CA.on.r.conf.level(n = n, n.median = n.median, 
            r = r, lpl.rank = lpl.rank, upl.rank = n + 1 - n.plus.one.minus.upl.rank, 
            pi.type = pi.type, integrate.args.list = integrate.args.list)
    })
    conf.level
}
