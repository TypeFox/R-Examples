ezmnorm <-
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
    n <- length(x)
    if (n < 1) 
        stop("'x' must contain at least one non-missing value")
    r <- sum(x == 0)
    if (r == 0) 
        warning("No 0 values in 'x'.")
    n <- length(x)
    phat <- r/n
    method <- match.arg(method)
    switch(method, mvue = {
        mean.zmnorm <- mean(x)
        if (r == n) {
            sd.zmnorm <- 0
            mean.norm <- NA
            sd.norm <- NA
        } else {
            x.non.zero <- x[x != 0]
            mean.norm <- mean(x.non.zero)
            sd.norm <- sd(x.non.zero)
            if (r == (n - 1)) {
                sd.zmnorm <- sqrt(x.non.zero^2/n)
            } else {
                sd.zmnorm <- sqrt(((n - r - 1)/(n - 1)) * sd.norm^2 + 
                  (r/n) * ((n - r)/(n - 1)) * mean.norm^2)
            }
        }
    })
    dist.params <- c(mean.norm, sd.norm, phat, mean.zmnorm, sd.zmnorm)
    names(dist.params) <- c("mean", "sd", "p.zero", "mean.zmnorm", 
        "sd.zmnorm")
    ret.list <- list(distribution = "Zero-Modified Normal", sample.size = n, 
        parameters = dist.params, n.param.est = 3, method = method, 
        data.name = data.name, bad.obs = bad.obs)
    if (ci) {
        if (n < 3) {
            warning(paste("Cannot compute confidence interval. ", 
                "You must have at least three non-missing observations."))
        }
        else {
            ci.type <- match.arg(ci.type, c("two-sided", "lower", 
                "upper"))
            if (conf.level <= 0 || conf.level >= 1) 
                stop("The value of 'conf.level' must be between 0 and 1")
            ci.obj <- ci.normal.approx(mean.zmnorm, sd.zmnorm/sqrt(n), 
                n = n, df = n - 2, ci.type = ci.type, alpha = 1 - 
                  conf.level)
            ci.obj$parameter <- "mean.zmnorm"
            ret.list <- c(ret.list, list(interval = ci.obj))
        }
    }
    oldClass(ret.list) <- "estimate"
    ret.list
}
