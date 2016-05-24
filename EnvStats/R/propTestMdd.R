propTestMdd <-
function (n.or.n1, n2 = n.or.n1, p0.or.p2 = 0.5, alpha = 0.05, 
    power = 0.95, sample.type = "one.sample", alternative = "two.sided", 
    two.sided.direction = "greater", approx = TRUE, correct = sample.type == 
        "two.sample", warn = TRUE, return.exact.list = TRUE, 
    tol = 1e-07, maxiter = 1000) 
{
    sample.type <- match.arg(sample.type, c("one.sample", "two.sample"))
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    two.sided.direction <- match.arg(two.sided.direction, c("greater", 
        "less"))
    if (!is.vector(n.or.n1, mode = "numeric") || !is.vector(p0.or.p2, 
        mode = "numeric") || !is.vector(alpha, mode = "numeric") || 
        !is.vector(power, mode = "numeric")) 
        stop("'n.or.n1', 'p0.or.p2', 'alpha', and 'power' must be numeric vectors.")
    if (!all(is.finite(n.or.n1)) || !all(is.finite(p0.or.p2)) || 
        !all(is.finite(alpha)) || !all(is.finite(power))) 
        stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
            "Undefined (Nan) values are not allowed in", "'n.or.n1', 'p0.or.p2', 'alpha', or 'power'"))
    if (any(n.or.n1 < 2)) 
        stop("All values of 'n.or.n1' must be greater than or equal to 2")
    if (!approx && !all(n.or.n1 == trunc(n.or.n1))) 
        stop("When approx=FALSE, all values of 'n.or.n1' must be integers")
    if (any(p0.or.p2 <= 0) || any(p0.or.p2 >= 1)) 
        stop("All values of 'p0.or.p2' must be greater than 0 and less than 1")
    if (any(alpha <= 0) || any(alpha >= 1)) 
        stop("All values of 'alpha' must be greater than 0 and less than 1")
    if (any(power < alpha) || any(power >= 1)) 
        stop(paste("All values of 'power' must be greater than or equal to", 
            "the corresponding elements of 'alpha', and less than 1"))
    if (sample.type == "two.sample" && !missing(n2)) {
        if (!is.vector(n2, mode = "numeric")) 
            stop("'n2' must be a numeric vector")
        if (!all(is.finite(n2))) 
            stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
                "Undefined (Nan) values are not allowed in 'n2'"))
        if (any(n2 < 2)) 
            stop("All values of 'n2' must be greater than or equal to 2")
        if (!approx) 
            stop(paste("Minimal detectable difference based on exact method", 
                "not available for two-sample test"))
    }
    arg.mat <- cbind.no.warn(n.or.n1 = as.vector(n.or.n1), n2 = as.vector(n2), 
        p0.or.p2 = as.vector(p0.or.p2), alpha = as.vector(alpha), 
        power = as.vector(power))
    for (i in c("n.or.n1", "n2", "p0.or.p2", "alpha", "power")) assign(i, 
        arg.mat[, i])
    N <- nrow(arg.mat)
    p.or.p1 <- numeric(N)
    index.0 <- power == alpha
    p.or.p1[index.0] <- p0.or.p2[index.0]
    extreme <- switch(alternative, less = .Machine$double.eps, 
        greater = 1 - .Machine$double.eps, two.sided = ifelse(two.sided.direction == 
            "less", .Machine$double.eps, 1 - .Machine$double.eps))
    extreme <- rep(extreme, N)
    power.extreme <- propTestPower(n.or.n1 = n.or.n1, p.or.p1 = extreme, 
        n2 = n2, p0.or.p2 = p0.or.p2, alpha = alpha, sample.type = sample.type, 
        alternative = alternative, approx = approx, correct = correct, 
        warn = FALSE, return.exact.list = FALSE)
    index.Inf <- power.extreme < power | is.na(power.extreme)
    if (any(index.Inf)) {
        p.or.p1[index.Inf] <- NA
        string <- ifelse(sample.type == "one.sample", "'n', 'p0', 'alpha', and 'power'", 
            "'n1', 'n2', 'p2', 'alpha', and 'power'")
        warning(paste("Elements with a missing value (NA) indicate", 
            "no attainable minimal detectable difference for the", 
            "given values of", string))
    }
    index.equal <- power.extreme == power
    p.or.p1[index.equal] <- extreme[index.equal]
    if (extreme[1] == .Machine$double.eps) {
        lower <- extreme
        upper <- p0.or.p2
        f.lower <- power - power.extreme
        f.upper <- power - alpha
    }
    else {
        lower <- p0.or.p2
        upper <- extreme
        f.lower <- power - alpha
        f.upper <- power - power.extreme
    }
    fcn.for.root <- function(p.or.p1, n.or.n1, n2, p0.or.p2, 
        alpha, power, sample.type, alternative, approx, correct) {
        power - propTestPower(n.or.n1 = n.or.n1, p.or.p1 = p.or.p1, 
            n2 = n2, p0.or.p2 = p0.or.p2, alpha = alpha, sample.type = sample.type, 
            alternative = alternative, approx = approx, correct = correct, 
            warn = FALSE, return.exact.list = FALSE)
    }
    for (i in (1:N)[!index.0 & !index.Inf & !index.equal]) {
        lower.i <- lower[i]
        upper.i <- upper[i]
        f.lower.i <- f.lower[i]
        f.upper.i <- f.upper[i]
        n.or.n1.i <- n.or.n1[i]
        p0.or.p2.i <- p0.or.p2[i]
        alpha.i <- alpha[i]
        power.i <- power[i]
        n2.i <- n2[i]
        p.or.p1[i] <- uniroot(fcn.for.root, lower = lower[i], 
            upper = upper.i, tol = tol, maxiter = maxiter, f.lower = f.lower.i, 
            f.upper = f.upper.i, n.or.n1 = n.or.n1.i, n2 = n2.i, 
            p0.or.p2 = p0.or.p2.i, alpha = alpha.i, power = power.i, 
            sample.type = sample.type, alternative = alternative, 
            approx = approx, correct = correct)$root
    }
    ret.val <- p.or.p1 - p0.or.p2
    names(ret.val) <- NULL
    if (approx && warn) {
        na.index <- is.na(p.or.p1)
        if (!all(na.index)) {
            if (sample.type == "one.sample" && any(index <- pmin(n.or.n1[!na.index] * 
                p.or.p1[!na.index], n.or.n1[!na.index] * (1 - 
                p.or.p1)[!na.index]) < 5)) 
                warning(paste("The sample size 'n' is too small,", 
                  "relative to the computed value of 'p', for the normal", 
                  "approximation to work well for the following", 
                  "element indices:\n\t", paste((1:N)[!na.index][index], 
                    collapse = " "), "\n\t"))
            if (sample.type == "two.sample" && any(index <- pmin(n.or.n1[!na.index] * 
                p.or.p1[!na.index], n.or.n1[!na.index] * (1 - 
                p.or.p1)[!na.index], n2[!na.index] * p0.or.p2[!na.index], 
                n2[!na.index] * (1 - p0.or.p2)[!na.index]) < 
                5)) 
                warning(paste("The sample sizes 'n1' and 'n2' are", 
                  "too small, relative to the computed value of", 
                  "'p1' and the given value of 'p2', for the normal", 
                  "approximation to work well for the following", 
                  "element indices:\n\t", paste((1:N)[!na.index][index], 
                    collapse = " "), "\n\t"))
        }
    }
    if (!approx && return.exact.list) {
        N <- length(p.or.p1)
        switch(alternative, less = {
            N.vars <- 3
            list.names <- c("power", "alpha", "q.critical.lower")
        }, greater = {
            N.vars <- 3
            list.names <- c("power", "alpha", "q.critical.upper")
        }, two.sided = {
            N.vars <- 4
            list.names <- c("power", "alpha", "q.critical.lower", 
                "q.critical.upper")
        })
        dum.list <- lapply(rep(N, N.vars), function(i) {
            x <- rep(NA, i)
            mode(x) <- "numeric"
            x
        })
        names(dum.list) <- list.names
        na.index <- is.na(p.or.p1)
        if (!all(na.index)) {
            power.list <- propTestPower(n.or.n1 = n.or.n1[!na.index], 
                p.or.p1 = p.or.p1[!na.index], n2 = n2[!na.index], 
                p0.or.p2 = p0.or.p2[!na.index], alpha = alpha[!na.index], 
                sample.type = sample.type, alternative = alternative, 
                approx = FALSE, warn = FALSE, return.exact.list = TRUE)
            dum.list <- lapply(1:N.vars, function(i, index, main.list, 
                sub.list) {
                main.list[[i]][index] <- sub.list[[i]]
                main.list[[i]]
            }, index = !na.index, main.list = dum.list, sub.list = power.list)
            names(dum.list) <- list.names
        }
        ret.val <- c(list(delta = ret.val), dum.list)
    }
    ret.val
}
