predIntNormN <-
function (half.width, n.mean = 1, k = 1, sigma.hat = 1, method = "Bonferroni", 
    conf.level = 0.95, round.up = TRUE, n.max = 5000, tol = 1e-07, 
    maxiter = 1000) 
{
    method <- match.arg(method, c("Bonferroni", "exact"))
    if (!is.vector(half.width, mode = "numeric") || !all(is.finite(half.width)) || 
        any(half.width <= .Machine$double.eps)) 
        stop(paste("'half.width' must be a numeric vector", "with all positive elements", 
            "and no Missing (NA), Infinite(-Inf, Inf),", "or Undefined (Nan) values."))
    if (!is.vector(k, mode = "numeric") || !all(is.finite(k)) || 
        !all(k == trunc(k)) || any(k < 1)) 
        stop(paste("'k' must be a numeric vector of integers", 
            "with all elements greater than 0", "and no Missing (NA), Infinite(-Inf, Inf),", 
            "or Undefined (Nan) values."))
    if (!is.vector(n.mean, mode = "numeric") || !all(is.finite(n.mean)) || 
        !all(n.mean == trunc(n.mean)) || any(n.mean < 1)) 
        stop(paste("'n.mean' must be a numeric vector of integers", 
            "with all elements greater than 0", "and no Missing (NA), Infinite(-Inf, Inf),", 
            "or Undefined (Nan) values."))
    if (!is.vector(sigma.hat, mode = "numeric") || !all(is.finite(sigma.hat)) || 
        any(sigma.hat < .Machine$double.eps)) 
        stop(paste("'sigma.hat' must be a numeric vector", "with all positive elements", 
            "and no Missing (NA), Infinite(-Inf, Inf),", "or Undefined (Nan) values."))
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
    arg.mat <- cbind.no.warn(half.width = as.vector(half.width), 
        n.mean = as.vector(n.mean), k = as.vector(k), sigma.hat = as.vector(sigma.hat), 
        conf.level = as.vector(conf.level))
    for (i in c("half.width", "n.mean", "k", "sigma.hat", "conf.level")) assign(i, 
        arg.mat[, i])
    N <- length(half.width)
    n.vec <- numeric(N)
    fcn.for.root <- function(n, hw, n.mean, k.weird, sigma.hat, 
        method, conf.level) {
        hw - predIntNormHalfWidth(n = n, n.mean = n.mean, k = k.weird, 
            sigma.hat = sigma.hat, method = method, conf.level = conf.level)
    }
    hw.2 <- predIntNormHalfWidth(n = 2, n.mean = n.mean, k = k, 
        sigma.hat = sigma.hat, method = method, conf.level = conf.level)
    hw.n.max <- predIntNormHalfWidth(n = n.max, n.mean = n.mean, 
        k = k, sigma.hat = sigma.hat, method = method, conf.level = conf.level)
    for (i in 1:N) {
        hw.i <- half.width[i]
        hw.2.i <- hw.2[i]
        if (hw.2.i <= hw.i) 
            n.vec[i] <- 2
        else {
            hw.n.max.i <- hw.n.max[i]
            if (hw.n.max.i > hw.i) {
                n.vec[i] <- NA
                warning(paste("Value of 'half.width' is too small for element ", 
                  i, ".  Try increasing the value of 'n.max'", 
                  sep = ""))
            }
            else {
                k.i <- k[i]
                n.mean.i <- n.mean[i]
                sigma.hat.i <- sigma.hat[i]
                conf.level.i <- conf.level[i]
                n.vec[i] <- uniroot(fcn.for.root, lower = 2, 
                  upper = n.max, f.lower = hw.i - hw.2.i, f.upper = hw.i - 
                    hw.n.max.i, hw = hw.i, n.mean = n.mean.i, 
                  k.weird = k.i, sigma.hat = sigma.hat.i, method = method, 
                  conf.level = conf.level.i, tol = tol, maxiter = maxiter)$root
            }
        }
    }
    if (round.up) 
        n.vec <- ceiling(n.vec)
    ret.val <- n.vec
    names(ret.val) <- NULL
    ret.val
}
