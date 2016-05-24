elogis <-
function (x, method = "mle", ci = FALSE, ci.type = "two-sided", 
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
    if (n < 2 || length(unique(x)) < 2) 
        stop(paste("'x' must contain at least 2 non-missing distinct values. ", 
            "This is not true for 'x' =", data.name))
    method <- match.arg(method, c("mle", "mme", "mmue"))
    location <- mean(x)
    scale <- (sqrt((n - 1)/n) * sd(x) * sqrt(3))/pi
    switch(method, mme = {
        dist.params <- c(location = location, scale = scale)
    }, mmue = {
        scale <- (sd(x) * sqrt(3))/pi
        dist.params <- c(location = location, scale = scale)
    }, mle = {
        neg.ll <- function(theta, y) {
            a <- theta[1]
            b <- theta[2]
            c <- (y - a)/b
            sum(c + log(b) + 2 * log(1 + exp(-c)))
        }
        dist.params <- nlminb(start = c(location, scale), objective = neg.ll, 
            lower = c(-Inf, .Machine$double.eps), y = x)$par
        names(dist.params) <- c("location", "scale")
    })
    ret.list <- list(distribution = "Logistic", sample.size = n, 
        parameters = dist.params, n.param.est = 2, method = method, 
        data.name = data.name, bad.obs = bad.obs)
    if (ci) {
        ci.type <- match.arg(ci.type, c("two-sided", "lower", 
            "upper"))
        ci.method <- match.arg(ci.method)
        if (conf.level <= 0 || conf.level >= 1) 
            stop("The value of 'conf.level' must be between 0 and 1.")
        ci.obj <- switch(ci.method, normal.approx = ci.normal.approx(theta.hat = dist.params["location"], 
            sd.theta.hat = (pi * dist.params["scale"])/sqrt(3 * 
                n), n, df = n - 1, ci.type = ci.type, alpha = 1 - 
                conf.level))
        ci.obj$parameter <- "location"
        ret.list <- c(ret.list, list(interval = ci.obj))
    }
    oldClass(ret.list) <- "estimate"
    ret.list
}
