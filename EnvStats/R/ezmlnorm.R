ezmlnorm <-
function (x, method = "mvue", ci = FALSE, ci.type = "two-sided", 
    ci.method = "normal.approx", conf.level = 0.95) 
{
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    data.name <- deparse(substitute(x))
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    if (any(x < 0)) 
        stop("All values of 'x' must be non-negative")
    r <- sum(x == 0)
    if (r == 0) 
        warning("No 0 values in 'x'.")
    n <- length(x)
    phat <- r/n
    method <- match.arg(method, "mvue")
    switch(method, mvue = {
        if (r == n) {
            mean.zmlnorm <- 0
            sd.zmlnorm <- 0
            mean.ln <- NA
            cv.ln <- NA
            mean.y <- NA
            sd.y <- NA
        } else {
            x.pos <- x[x > 0]
            y <- log(x.pos)
            mean.y <- mean(y)
            sd.y <- sd(y)
            if (r == (n - 1)) {
                mean.zmlnorm <- x.pos/n
                sd.zmlnorm <- x.pos/sqrt(n)
                mean.ln <- x.pos
                cv.ln <- NA
            } else {
                params.ln <- elnormAlt(x.pos)$parameters
                mean.ln <- params.ln["mean"]
                cv.ln <- params.ln["cv"]
                mean.zmlnorm <- (1 - phat) * mean.ln
                s2 <- sd.y^2
                sd.zmlnorm <- sqrt((1 - phat) * exp(2 * mean.y) * 
                  (finneys.g(n - r - 1, 2 * s2) - ((1 - r/(n - 
                    1)) * finneys.g(n - r - 1, ((n - r - 2) * 
                    s2)/(n - r - 1)))))
            }
        }
    })
    dist.params <- c(mean.y, sd.y, phat, mean.zmlnorm, sd.zmlnorm)
    names(dist.params) <- c("meanlog", "sdlog", "p.zero", "mean.zmlnorm", 
        "sd.zmlnorm")
    ret.list <- list(distribution = "Zero-Modified Lognormal (Delta)", 
        sample.size = n, parameters = dist.params, n.param.est = 3, 
        method = method, data.name = data.name, bad.obs = bad.obs)
    if (ci) {
        if (n < 3 || r == n) 
            warning(paste("Cannot compute confidence interval. ", 
                "You must have at least three non-missing observations,", 
                "and at least one observation must be non-zero."))
        else {
            ci.method <- match.arg(ci.method, "normal.approx")
            ci.type <- match.arg(ci.type, c("two-sided", "lower", 
                "upper"))
            if (conf.level <= 0 || conf.level >= 1) 
                stop("The value of 'conf.level' must be between 0 and 1")
            avar.mean.zmlnorm <- (exp(2 * mean.y + sd.y^2) * 
                (1 - phat) * (phat + 0.5 * (2 * sd.y^2 + sd.y^4)))/n
            ci.obj <- ci.normal.approx(mean.zmlnorm, sqrt(avar.mean.zmlnorm), 
                n = n, df = n - 2, ci.type = ci.type, alpha = 1 - 
                  conf.level, lb = 0)
            ci.obj$parameter <- "mean.zmlnorm"
            ret.list <- c(ret.list, list(interval = ci.obj))
        }
    }
    oldClass(ret.list) <- "estimate"
    ret.list
}
