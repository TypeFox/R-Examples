eevd <-
function (x, method = "mle", pwme.method = "unbiased", plot.pos.cons = c(a = 0.35, 
    b = 0), ci = FALSE, ci.parameter = "location", ci.type = "two-sided", 
    ci.method = "normal.approx", conf.level = 0.95) 
{
    if (!is.vector(x, mode = "numeric")) 
        stop("'x' must be a numeric vector")
    data.name <- deparse(substitute(x))
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    n <- length(x)
    if (n < 2 || length(unique(x)) < 2) 
        stop("'x' must contain at least 2 non-missing distinct values")
    method <- match.arg(method, c("mle", "mme", "mmue", "pwme"))
    ret.list.method <- method
    ci.parameter <- match.arg(ci.parameter, c("location", "scale"))
    mean.x <- mean(x)
    sd.x <- sd(x)
    scale.mmue <- (sd.x * sqrt(6))/pi
    location.mmue <- mean.x - .Eulers.constant * scale.mmue
    switch(method, mmue = {
        dist.params <- c(location = location.mmue, scale = scale.mmue)
    }, mme = {
        scale.mme <- sqrt((n - 1)/n) * scale.mmue
        location.mme <- mean.x - .Eulers.constant * scale.mme
        dist.params <- c(location = location.mme, scale = scale.mme)
    }, mle = {
        fcn.to.min <- function(scale, x, x.bar) {
            vec <- exp(-x/scale)
            (scale - (x.bar - sum(x * vec)/sum(vec)))^2
        }
        scale.mle <- nlminb(start = scale.mmue, objective = fcn.to.min, 
            lower = .Machine$double.eps, x = x, x.bar = mean.x)$par
        location.mle <- -scale.mle * log(mean(exp(-x/scale.mle)))
        dist.params <- c(location = location.mle, scale = scale.mle)
    }, pwme = {
        pwme.method <- match.arg(pwme.method, c("unbiased", "plotting.position"))
        M0 <- mean.x
        scale.pwme <- (M0 - 2 * pwMoment(x, k = 1, method = pwme.method, 
            plot.pos.cons = plot.pos.cons))/log(2)
        location.pwme <- M0 - .Eulers.constant * scale.pwme
        dist.params <- c(location = location.pwme, scale = scale.pwme)
        if (pwme.method == "unbiased") ret.list.method <- paste("Unbiased", 
            method) else ret.list.method <- paste("Plotting-position ", 
            method, "\n", space(33), "with coefficients\n", space(33), 
            "(a, b) = (", plot.pos.cons["a"], ", ", plot.pos.cons["b"], 
            ")", sep = "")
    })
    ret.list <- list(distribution = "Extreme Value", sample.size = n, 
        parameters = dist.params, n.param.est = 2, method = ret.list.method, 
        data.name = data.name, bad.obs = bad.obs)
    if (ci) {
        ci.type <- match.arg(ci.type, c("two-sided", "lower", 
            "upper"))
        ci.method <- match.arg(ci.method)
        if (conf.level <= 0 || conf.level >= 1) 
            stop("The value of 'conf.level' must be between 0 and 1.")
        scale <- dist.params["scale"]
        if (method == "pwme") {
            L <- dilog(1.5)
            pso6 <- pi^2/6
            pso4 <- pi^2/4
            l2 <- log(2)
            l22 <- l2^2
            Var.v <- pso6 * n
            Var.w <- (pso6 - l22 - 2 * L) * n^3 + (-pso4 + l22 + 
                6 * L) * n^2 + (pso4 - 4 * L) * n
            Cov.vw <- (pso6/2 - l22/2) * n^2 + (pso6/2 + l22/2) * 
                n
            den <- n * (n - 1) * l2
        }
        sd.param <- switch(ci.parameter, location = switch(method, 
            mmue = , mme = sqrt((1.167814 * scale^2)/n), mle = sqrt((1.10867 * 
                scale^2)/n), pwme = {
                c1 <- 1/n
                c2 <- ((n + 1) * .Eulers.constant)/den
                c3 <- (2 * .Eulers.constant)/den
                scale * sqrt((c1 + c2)^2 * Var.v + c3^2 * Var.w - 
                  2 * (c1 + c2) * c3 * Cov.vw)
            }), scale = switch(method, mmue = , mme = sqrt((1.1 * 
            scale^2)/n), mle = sqrt(((6/pi^2) * scale^2)/n), 
            pwme = {
                c1 <- (n + 1)/den
                c2 <- 2/den
                scale * sqrt(c1^2 * Var.v + c2^2 * Var.w - 2 * 
                  c1 * c2 * Cov.vw)
            }))
        ci.obj <- switch(ci.method, normal.approx = ci.normal.approx(dist.params[ci.parameter], 
            sd.param, n = n, df = n - 1, ci.type = ci.type, alpha = 1 - 
                conf.level, lb = ifelse(ci.parameter == "scale", 
                0, -Inf)))
        ci.obj$parameter <- ci.parameter
        ret.list <- c(ret.list, list(interval = ci.obj))
    }
    oldClass(ret.list) <- "estimate"
    ret.list
}
