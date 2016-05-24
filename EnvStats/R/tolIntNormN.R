tolIntNormN <-
function (half.width, sigma.hat = 1, coverage = 0.95, cov.type = "content", 
    conf.level = 0.95, method = "wald.wolfowitz", round.up = TRUE, 
    n.max = 5000, tol = 1e-07, maxiter = 1000) 
{
    cov.type <- match.arg(cov.type, c("content", "expectation"))
    method <- match.arg(method, c("exact", "wald.wolfowitz"))
    if (method != "wald.wolfowitz") 
        stop("Currently only method='wald.wolfowitz' is allowed.")
    if (!is.vector(half.width, mode = "numeric") || !all(is.finite(half.width)) || 
        any(half.width <= .Machine$double.eps)) 
        stop(paste("'half.width' must be a numeric vector", "with all positive elements", 
            "and no Missing (NA), Infinite(-Inf, Inf),", "or Undefined (Nan) values."))
    if (!is.vector(sigma.hat, mode = "numeric") || !all(is.finite(sigma.hat)) || 
        any(sigma.hat < .Machine$double.eps)) 
        stop(paste("'sigma.hat' must be a numeric vector", "with all positive elements", 
            "and no Missing (NA), Infinite(-Inf, Inf),", "or Undefined (Nan) values."))
    if (!is.vector(coverage, mode = "numeric") || !all(is.finite(coverage)) || 
        any(coverage <= .Machine$double.eps) || any(coverage >= 
        1 - .Machine$double.eps)) 
        stop(paste("'coverage' must be a numeric vector", "with all elements between 0 and 1", 
            "and no Missing (NA), Infinite(-Inf, Inf),", "or Undefined (Nan) values."))
    if (!is.vector(n.max, mode = "numeric") || length(n.max) != 
        1 || !is.finite(n.max) || n.max != trunc(n.max) || n.max < 
        2) 
        stop("'n.max' must be a positive integer greater than 1")
    if (!is.vector(maxiter, mode = "numeric") || length(maxiter) != 
        1 || !is.finite(maxiter) || maxiter != trunc(maxiter) || 
        maxiter < 2) 
        stop("'maxiter' must be a positive integer greater than 1")
    if (cov.type == "content") {
        if (!is.vector(conf.level, mode = "numeric") || !all(is.finite(conf.level)) || 
            any(conf.level <= .Machine$double.eps) || any(conf.level >= 
            1 - .Machine$double.eps)) 
            stop(paste("'conf.level' must be a numeric vector", 
                "with all elements between 0 and 1", "and no Missing (NA), Infinite(-Inf, Inf),", 
                "or Undefined (Nan) values."))
        arg.mat <- cbind.no.warn(half.width = as.vector(half.width), 
            sigma.hat = as.vector(sigma.hat), coverage = as.vector(coverage), 
            conf.level = as.vector(conf.level))
        for (i in c("half.width", "sigma.hat", "coverage", "conf.level")) assign(i, 
            arg.mat[, i])
        N <- length(half.width)
        n.vec <- numeric(N)
        fcn.for.root <- function(n, hw, sigma.hat, coverage, 
            conf.level, method) {
            hw - tolIntNormHalfWidth(n = n, sigma.hat = sigma.hat, 
                coverage = coverage, cov.type = "content", conf.level = conf.level, 
                method = method)
        }
        hw.2 <- tolIntNormHalfWidth(n = 2, sigma.hat = sigma.hat, 
            coverage = coverage, cov.type = "content", conf.level = conf.level, 
            method = method)
        hw.n.max <- tolIntNormHalfWidth(n = n.max, sigma.hat = sigma.hat, 
            coverage = coverage, cov.type = "content", conf.level = conf.level, 
            method = method)
        for (i in 1:N) {
            hw.i <- half.width[i]
            hw.2.i <- hw.2[i]
            hw.n.max.i <- hw.n.max[i]
            if (hw.2.i <= hw.i) {
                n.vec[i] <- 2
            }
            else {
                if (hw.n.max.i > hw.i) {
                  n.vec[i] <- NA
                  warning(paste("Value of 'half.width' is too small", 
                    "for element", i, ".  Try increasing the value of 'n.max'.", 
                    sep = ""))
                }
                else {
                  sigma.hat.i <- sigma.hat[i]
                  coverage.i <- coverage[i]
                  conf.level.i <- conf.level[i]
                  n.vec[i] <- uniroot(fcn.for.root, lower = 2, 
                    upper = n.max, f.lower = hw.i - hw.2.i, f.upper = hw.i - 
                      hw.n.max.i, hw = hw.i, sigma.hat = sigma.hat.i, 
                    coverage = coverage.i, conf.level = conf.level.i, 
                    method = method, tol = tol, maxiter = maxiter)$root
                }
            }
        }
        if (round.up) 
            n.vec <- ceiling(n.vec)
        ret.val <- n.vec
        names(ret.val) <- NULL
    }
    else {
        ret.val <- predIntNormN(half.width = half.width, k = 1, 
            n.mean = 1, sigma.hat = sigma.hat, conf.level = coverage, 
            round.up = round.up, n.max = n.max, tol = tol, maxiter = maxiter)
    }
    ret.val
}
