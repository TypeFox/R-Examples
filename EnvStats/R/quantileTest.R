quantileTest <-
function (x, y, alternative = "greater", target.quantile = 0.5, 
    target.r = NULL, exact.p = TRUE) 
{
    if (!is.vector(x, mode = "numeric") || is.factor(x) || !is.vector(y, 
        mode = "numeric") || is.factor(y)) 
        stop("'x' and 'y' must be numeric vectors")
    data.name <- c(deparse(substitute(x)), deparse(substitute(y)))
    names(data.name) <- c("x", "y")
    if ((bad.obs.x <- sum(!(x.ok <- is.finite(x)))) > 0) {
        x <- x[x.ok]
        warning(paste(bad.obs.x, "observations with NA/NaN/Inf in 'x' removed."))
    }
    nx <- length(x)
    if (nx < 2) 
        stop(paste("'x' must contain at least two non-missing,", 
            "finite observations"))
    if ((bad.obs.y <- sum(!(y.ok <- is.finite(y)))) > 0) {
        is.not.finite.warning(y)
        y <- y[y.ok]
        warning(paste(bad.obs.y, "observations with NA/NaN/Inf in 'y' removed."))
    }
    ny <- length(y)
    if (ny < 2) 
        stop(paste("'y' must contain at least two non-missing,", 
            "finite observations"))
    alternative <- match.arg(alternative, c("greater", "less"))
    N <- nx + ny
    vec <- c(x, y)
    index <- c(rep(1, nx), rep(2, ny))
    sorted.vec <- sort(vec)
    sorted.index <- index[order(vec)]
    sorted.rank.vec <- rank(sorted.vec)
    if (is.null(target.r)) {
        if (!is.vector(target.quantile, mode = "numeric") || 
            length(target.quantile) != 1 || target.quantile <= 
            0 || target.quantile >= 1) 
            stop(paste("'target.quantile' must be a numeric scalar", 
                "between 0 and 1"))
        b <- target.quantile
        threshold.rank <- (N + 1) * b
        ceiling.tr <- ceiling(threshold.rank)
        threshold.rank <- ifelse(ceiling.tr - threshold.rank < 
            .Machine$double.eps, ceiling.tr + 1, ceiling.tr)
        if (threshold.rank <= 1) 
            stop(paste("The value of 'target.quantile' is", "too small relative to the sample size.  You must", 
                "increase the value of 'target.quantile' to at", 
                "least", 2/(N + 1)))
        if (threshold.rank > N) 
            stop(paste("The value of 'target.quantile' is", "too large relative to the sample size.  You must", 
                "decrease the value of 'target.quantile' to at", 
                "most", N/(N + 1)))
    }
    else {
        if (!is.vector(target.r, mode = "numeric") || length(target.r) != 
            1 || target.r != trunc(target.r) || target.r <= 1 || 
            target.r > N) 
            stop(paste("'target.r' must be a positive integer strictly", 
                "bigger than 1 and less than or equal to", N))
        threshold.rank <- sorted.rank.vec[N - target.r + 1]
    }
    quantile.ub <- threshold.rank/(N + 1)
    index <- sorted.rank.vec >= threshold.rank
    r <- sum(index)
    if (alternative == "greater") {
        k <- sum(index & sorted.index == 1)
        m <- nx
        n <- ny
    }
    else {
        k <- sum(index & sorted.index == 2)
        m <- ny
        n <- nx
    }
    if (exact.p) {
        p.value <- 1 - phyper(q = k - 1, m = m, n = n, k = r)
        method <- "Quantile Test"
    }
    else {
        p.value <- 1 - pnorm((k - (m * r)/N - 0.5)/sqrt((nx * 
            ny * r * (N - r))/(N^2 * (N - 1))))
        method <- "Quantile Test (approximate p-value)"
    }
    stat <- c(k, r)
    if (alternative == "greater") {
        names(stat) <- c("k (# x obs of r largest)", "r")
    }
    else {
        names(stat) <- c("k (# y obs of r largest)", "r")
    }
    parameters <- c(m = nx, n = ny, quantile.ub = quantile.ub)
    null.value <- 0
    names(null.value) <- "e"
    if (alternative == "greater") {
        alternative <- paste("Tail of Fx Shifted to Right of", 
            "Tail of Fy.", "0 < e <= 1, where", "Fx(t) = (1-e)*Fy(t) + e*Fz(t),", 
            "Fz(t) <= Fy(t) for all t,", "and Fy != Fz", sep = paste("\n", 
                space(33), sep = ""))
    }
    else {
        alternative <- paste("Tail of Fx Shifted to Left of", 
            "Tail of Fy.", "0 < e <= 1, where", "Fy(t) = (1-e)*Fx(t) + e*Fz(t),", 
            "Fz(t) <= Fx(t) for all t,", "and Fx != Fz", sep = paste("\n", 
                space(33), sep = ""))
    }
    ret.list <- list(statistic = stat, parameters = parameters, 
        p.value = p.value, estimate = NULL, null.value = null.value, 
        alternative = alternative, method = method, estimation.method = NULL, 
        sample.size = c(nx = nx, ny = ny), data.name = data.name, 
        bad.obs = bad.obs.x + bad.obs.y)
    oldClass(ret.list) <- "htest"
    ret.list
}
