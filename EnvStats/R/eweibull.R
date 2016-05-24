eweibull <-
function (x, method = "mle") 
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
            "and all non-missing values of 'x' must be non-negative. ", 
            "This is not true for 'x' =", data.name))
    method <- match.arg(method, c("mle", "mme", "mmue"))
    x.bar <- mean(x)
    sd.x <- ifelse(method == "mmue", sd(x), sqrt((n - 1)/n) * 
        sd(x))
    mcf <- function(c, x.bar, sd.x) {
        t1 <- gamma((c + 2)/c)
        t2 <- gamma((c + 1)/c)
        t3 <- sqrt((t1/(t2^2)) - 1)
        ((sd.x/x.bar) - t3)^2
    }
    shape <- nlminb(start = 1, objective = mcf, lower = .Machine$double.eps, 
        x.bar = x.bar, sd.x = sd.x)$par
    scale <- x.bar/gamma((shape + 1)/shape)
    switch(method, mme = , mmue = {
        dist.params <- c(shape = shape, scale = scale)
    }, mle = {
        neg.ll <- function(theta, y) {
            b <- theta[1]
            c <- theta[2]
            sum(-log(c) - (c - 1) * log(y) + c * log(b) + ((y/b)^c))
        }
        theta.hat <- nlminb(start = c(scale, shape), objective = neg.ll, 
            lower = .Machine$double.eps, y = x)$par
        dist.params <- c(shape = theta.hat[2], scale = theta.hat[1])
    })
    ret.list <- list(distribution = "Weibull", sample.size = n, 
        parameters = dist.params, n.param.est = 2, method = method, 
        data.name = data.name, bad.obs = bad.obs)
    oldClass(ret.list) <- "estimate"
    ret.list
}
