ciNormN <-
function (half.width, sigma.hat = 1, conf.level = 0.95, sample.type = ifelse(is.null(n2), 
    "one.sample", "two.sample"), n2 = NULL, round.up = TRUE, 
    n.max = 5000, tol = 1e-07, maxiter = 1000) 
{
    sample.type <- match.arg(sample.type, c("one.sample", "two.sample"))
    if (!is.vector(half.width, mode = "numeric") || !is.vector(sigma.hat, 
        mode = "numeric") || !is.vector(conf.level, mode = "numeric")) 
        stop("'half.width', 'sigma.hat', and 'conf.level' must be numeric vectors.")
    if (!all(is.finite(half.width)) || !all(is.finite(sigma.hat)) || 
        !all(is.finite(conf.level))) 
        stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
            "Undefined (Nan) values are not allowed in", "'half.width', 'sigma.hat', or 'conf.level'"))
    if (any(half.width < .Machine$double.eps)) 
        stop("All values of 'half.width' must be positive")
    if (any(sigma.hat < .Machine$double.eps)) 
        stop("All values of 'sigma.hat' must be positive.")
    if (any(conf.level <= .Machine$double.eps) || any(conf.level >= 
        1 - .Machine$double.eps)) 
        stop("All values of 'conf.level' must be greater than 0 and less than 1")
    if (!is.vector(n.max, mode = "numeric") || length(n.max) != 
        1 || !is.finite(n.max) || n.max != trunc(n.max) || n.max < 
        2) 
        stop("'n.max' must be a positive integer greater than 1")
    if (!is.vector(maxiter, mode = "numeric") || length(maxiter) != 
        1 || !is.finite(maxiter) || maxiter != trunc(maxiter) || 
        maxiter < 2) 
        stop("'maxiter' must be a positive integer greater than 1")
    if (n2.constrained <- sample.type == "two.sample" && !is.null(n2)) {
        if (!is.vector(n2, mode = "numeric")) 
            stop("'n2' must be a numeric vector")
        if (!all(is.finite(n2))) 
            stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
                "Undefined (Nan) values are not allowed in 'n2'"))
        if (any(n2 != trunc(n2)) || any(n2 < 2)) 
            stop("All values of 'n2' must be positive integers larger than 1")
        arg.mat <- cbind.no.warn(half.width = as.vector(half.width), 
            sigma.hat = as.vector(sigma.hat), conf.level = as.vector(conf.level), 
            n2 = as.vector(n2))
        for (i in c("half.width", "sigma.hat", "conf.level", 
            "n2")) assign(i, arg.mat[, i])
    }
    else {
        arg.mat <- cbind.no.warn(half.width = as.vector(half.width), 
            sigma.hat = as.vector(sigma.hat), conf.level = as.vector(conf.level))
        for (i in c("half.width", "sigma.hat", "conf.level")) assign(i, 
            arg.mat[, i])
    }
    N <- length(half.width)
    n.vec <- numeric(N)
    type.fac <- ifelse(sample.type == "two.sample", 2, 1)
    sigma.hat.fac <- sqrt(type.fac) * sigma.hat
    sohw2 <- (sigma.hat.fac/half.width)^2
    crit.p <- 1 - (1 - conf.level)/2
    fcn.for.root <- function(n, hw, sigma.hat.fac, crit.p, type.fac) {
        hw - (sigma.hat.fac/sqrt(n)) * qt(crit.p, type.fac * 
            (n - 1))
    }
    hw.2 <- (sigma.hat.fac/sqrt(2)) * qt(crit.p, type.fac)
    hw.n.max <- (sigma.hat.fac/sqrt(n.max)) * qt(crit.p, type.fac * 
        (n.max - 1))
    for (i in 1:N) {
        hw.i <- half.width[i]
        hw.2.i <- hw.2[i]
        if (hw.2.i <= hw.i) 
            n.vec[i] <- 2
        else {
            hw.n.max.i <- hw.n.max[i]
            if (hw.n.max.i > hw.i) {
                n.vec[i] <- NA
                warning(paste("Error in algorithm for element ", 
                  i, ".\n", "Try increasing the value of the argument 'n.max'", 
                  sep = ""))
            }
            else {
                sigma.hat.fac.i <- sigma.hat.fac[i]
                crit.p.i <- crit.p[i]
                n.vec[i] <- uniroot(fcn.for.root, lower = 2, 
                  upper = n.max, f.lower = hw.i - hw.2.i, f.upper = hw.i - 
                    hw.n.max.i, hw = hw.i, sigma.hat.fac = sigma.hat.fac.i, 
                  crit.p = crit.p.i, type.fac = type.fac, tol = tol, 
                  maxiter = maxiter)$root
            }
        }
    }
    if (n2.constrained) {
        n1 <- (n.vec * n2)/(2 * n2 - n.vec)
        if (any(index <- !is.finite(n1) | n1 < 0)) {
            n1[index] <- NA
            warning(paste("One or more constrained values of 'n2' is(are)", 
                "too small given the associated values of", "'half.width', 'sigma.hat', and 'conf.level'"))
        }
        names(n1) <- names(n2) <- NULL
        if (round.up) 
            n1 <- ceiling(n1)
        ret.val <- list(n1 = n1, n2 = n2)
    }
    else {
        if (round.up) 
            n.vec <- ceiling(n.vec)
        ret.val <- n.vec
        names(ret.val) <- NULL
    }
    ret.val
}
