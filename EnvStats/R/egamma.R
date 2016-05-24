egamma <-
function (x, method = "mle", ci = FALSE, ci.type = "two-sided", 
    ci.method = "normal.approx", normal.approx.transform = "kulkarni.powar", 
    conf.level = 0.95) 
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
    if (n < 2 || any(x < 0) || length(unique(x)) < 2) 
        stop(paste("'x' must contain at least 2 non-missing distinct values,", 
            "and all non-missing values of x must be non-negative."))
    method <- match.arg(method, c("mle", "bcmle", "mme", "mmue"))
    m <- mean(x)
    s <- sqrt((n - 1)/n) * sd(x)
    shape <- (m/s)^2
    if (method != "mme") {
        shape <- switch(method, mmue = (m/(sqrt(n/(n - 1)) * 
            s))^2, bcmle = , mle = {
            msf <- function(shape, lmx, mlx) {
                (log(shape) - digamma(shape) - lmx + mlx)^2
            }
            nlminb(start = shape, objective = msf, lower = .Machine$double.eps, 
                lmx = log(m), mlx = mean(log(x)))$par
        })
    }
    if (method == "bcmle") {
        shape <- ((n - 3)/n) * shape + (2/(3 * n))
    }
    scale <- m/shape
    method.string <- switch(method, mle = "MLE", bcmle = "Bias-Corrected MLE", 
        mme = "Method of Moments", mmue = paste("Method of Moments Based on\n", 
            space(33), "Unbiased Variance Estimate", sep = ""))
    ret.list <- list(distribution = "Gamma", sample.size = n, 
        parameters = c(shape = shape, scale = scale), method = method.string, 
        data.name = data.name, bad.obs = bad.obs)
    if (ci) {
        ci.type <- match.arg(ci.type, c("two-sided", "lower", 
            "upper"))
        ci.method <- match.arg(ci.method, c("normal.approx", 
            "chisq.approx", "profile.likelihood"))
        if (ci.method == "profile.likelihood") {
            if (method != "mle") 
                stop("When ci.method=\"profile.likelihood\" you must set method=\"mle\"")
            if (ci.type != "two-sided") 
                stop("When ci.method=\"profile.likelihood\" you must set ci.type=\"two-sided\"")
        }
        normal.approx.transform <- match.arg(normal.approx.transform, 
            c("kulkarni.powar", "cube.root", "fourth.root"))
        if (conf.level <= 0 || conf.level >= 1) 
            stop("The value of 'conf.level' must be between 0 and 1.")
        if (ci.method %in% c("normal.approx", "profile.likelihood")) {
            ci.obj <- ci.gamma.normal.approx(x = x, shape = shape, 
                shape.est.method = method, ci.type = ci.type, 
                conf.level = conf.level, normal.approx.transform = normal.approx.transform)
            if (ci.method == "profile.likelihood") {
                limits <- ci.obj$limits
                names(limits) <- NULL
                ci.obj <- ci.gamma.profile.likelihood(x = x, 
                  shape.mle = shape, scale.mle = scale, conf.level = conf.level, 
                  LCL.start = limits[1], UCL.start = limits[2])
            }
        }
        else {
            ci.obj <- ci.gamma.chisq.approx(x = x, shape = shape, 
                shape.est.method = method, ci.type = ci.type, 
                conf.level = conf.level)
        }
        ret.list <- c(ret.list, list(interval = ci.obj))
    }
    oldClass(ret.list) <- "estimate"
    ret.list
}
