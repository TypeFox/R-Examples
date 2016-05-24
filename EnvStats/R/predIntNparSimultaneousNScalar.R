predIntNparSimultaneousNScalar <-
function (n.median = 1, k = 1, m = 2, r = 1, rule = c("k.of.m", 
    "CA", "Modified.CA"), lpl.rank = ifelse(pi.type == "upper", 
    0, 1), n.plus.one.minus.upl.rank = ifelse(pi.type == "lower", 
    0, 1), pi.type = "upper", conf.level = 0.95, n.max = 5000, 
    integrate.args.list = NULL, maxiter = 1000) 
{
    if (!is.vector(n.median, mode = "numeric") || length(n.median) != 
        1 || n.median != trunc(n.median) || n.median < 1 || !is.odd(n.median)) 
        stop("'n.median' must be a positive odd integer")
    rule <- match.arg(rule)
    switch(rule, k.of.m = {
        if (!is.vector(k, mode = "numeric") || length(k) != 1 || 
            k != trunc(k) || k < 1 || !is.vector(m, mode = "numeric") || 
            length(m) != 1 || m != trunc(m) || m < 1 || !is.vector(r, 
            mode = "numeric") || length(r) != 1 || r != trunc(r) || 
            r < 1 || k > m) stop(paste("'k', 'm', and 'r' must be positive integers,", 
            "and 'k' must be between 1 and 'm'"))
        if (k > 1) stop("Finding n for the k-of-m rule with k>1 is not yet implemented")
    }, CA = {
        if (!is.vector(m, mode = "numeric") || length(m) != 1 || 
            m != trunc(m) || m < 1 || !is.vector(r, mode = "numeric") || 
            length(r) != 1 || r != trunc(r) || r < 1) stop("'m', and 'r' must be positive integers")
        k <- "First.or.all.of.next.m.minus.one"
    }, Modified.CA = {
        if (!is.vector(m, mode = "numeric") || length(m) != 1 || 
            m != trunc(m) || m < 1 || !is.vector(r, mode = "numeric") || 
            length(r) != 1 || r != trunc(r) || r < 1) stop("'m', and 'r' must be positive integers")
        k <- "First.or.at.least.two.of.next.three"
        m <- 4
    })
    if (!is.vector(n.max, mode = "numeric") || length(n.max) != 
        1 || !is.finite(n.max) || n.max != trunc(n.max) || n.max < 
        2) 
        stop("'n.max' must be a positive integer greater than 1")
    pi.type <- match.arg(pi.type, c("upper", "lower"))
    if (pi.type == "upper") 
        lpl.rank <- 0
    else n.plus.one.minus.upl.rank <- 0
    if (!is.vector(lpl.rank, mode = "numeric") || length(lpl.rank) > 
        1 || !is.finite(lpl.rank) || lpl.rank != trunc(lpl.rank) || 
        lpl.rank < 0 || lpl.rank >= n.max) 
        stop("'lpl.rank' must be a non-negative integer less than 'n.max'")
    if (pi.type == "lower" & lpl.rank < 1) 
        stop("When pi.type='lower', 'lpl.rank' must be a positive integer")
    if (!is.vector(n.plus.one.minus.upl.rank, mode = "numeric") || 
        length(n.plus.one.minus.upl.rank) > 1 || !is.finite(n.plus.one.minus.upl.rank) || 
        n.plus.one.minus.upl.rank != trunc(n.plus.one.minus.upl.rank) || 
        n.plus.one.minus.upl.rank < 0 || n.plus.one.minus.upl.rank >= 
        n.max) 
        stop(paste("'n.plus.one.minus.upl.rank' must be a non-negative", 
            "integer less than 'n.max'"))
    if (pi.type == "upper" & n.plus.one.minus.upl.rank < 1) 
        stop("When pi.type='upper', 'n.plus.one.minus.upl.rank' must be positive integer")
    min.n <- ifelse(pi.type == "lower", max(2, lpl.rank + 1), 
        max(2, n.plus.one.minus.upl.rank + 1))
    conf.level.min.n <- predIntNparSimultaneousConfLevelScalar(n = min.n, 
        n.median = n.median, k = k, m = m, r = r, rule = rule, 
        lpl.rank = lpl.rank, n.plus.one.minus.upl.rank = n.plus.one.minus.upl.rank, 
        pi.type = pi.type, integrate.args.list = integrate.args.list)
    if (conf.level.min.n >= conf.level) {
        n <- min.n
    }
    else if (rule == "k.of.m" && r == 1 && n.median == 1) {
        n <- predIntNparNScalar(k = k, m = m, lpl.rank = lpl.rank, 
            n.plus.one.minus.upl.rank = n.plus.one.minus.upl.rank, 
            pi.type = pi.type, conf.level = conf.level, n.max = n.max, 
            maxiter = maxiter)
    }
    else {
        conf.level.n.max <- predIntNparSimultaneousConfLevelScalar(n = n.max, 
            n.median = n.median, k = k, m = m, r = r, rule = rule, 
            lpl.rank = lpl.rank, n.plus.one.minus.upl.rank = n.plus.one.minus.upl.rank, 
            pi.type = pi.type, integrate.args.list = integrate.args.list)
        if (conf.level.n.max < conf.level) {
            stop(paste("The required confidence level cannot be achieved", 
                "with the supplied value of 'n.max'.  You need to", 
                "adjust the value of 'conf.level', 'n.median', 'k', 'm', 'r',", 
                "and/or 'n.max'"))
        }
        fcn.for.root <- function(n.weird, n.median.weird, k.weird, 
            m.weird, r.weird, rule.weird, lpl.rank.weird, n.plus.one.minus.upl.rank.weird, 
            pi.type.weird, conf.level.weird, integrate.args.list.weird) {
            n <- trunc(n.weird)
            predIntNparSimultaneousConfLevelScalar(n = n, n.median = n.median.weird, 
                k = k.weird, m = m.weird, r = r.weird, rule = rule.weird, 
                lpl.rank = lpl.rank.weird, n.plus.one.minus.upl.rank = n.plus.one.minus.upl.rank.weird, 
                pi.type = pi.type.weird, integrate.args.list = integrate.args.list.weird) - 
                conf.level.weird
        }
        uni.list <- uniroot(fcn.for.root, interval = c(min.n, 
            n.max), tol = 1, maxiter = maxiter, n.median.weird = n.median, 
            k.weird = k, m.weird = m, r.weird = r, rule.weird = rule, 
            lpl.rank.weird = lpl.rank, n.plus.one.minus.upl.rank.weird = n.plus.one.minus.upl.rank, 
            pi.type.weird = pi.type, conf.level.weird = conf.level, 
            integrate.args.list.weird = integrate.args.list)
        if (uni.list$f.root > 0) 
            n <- trunc(uni.list$root)
        else n <- trunc(uni.list$root) + 1
    }
    n
}
