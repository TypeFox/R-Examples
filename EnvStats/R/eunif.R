eunif <-
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
    if (n < 2 || length(unique(x)) < 2) 
        stop(paste("'x' must contain at least 2 non-missing distinct values. ", 
            "This is not true for 'x' =", data.name))
    method <- match.arg(method, c("mle", "mme", "mmue"))
    switch(method, mme = {
        a <- mean(x)
        h <- sqrt(3) * sqrt((n - 1)/n) * sd(x)
        dist.params <- c(min = a - h, max = a + h)
    }, mmue = {
        a <- mean(x)
        h <- sqrt(3) * sd(x)
        dist.params <- c(min = a - h, max = a + h)
    }, mle = {
        dist.params <- c(min = min(x), max = max(x))
    })
    ret.list <- list(distribution = "Uniform", sample.size = n, 
        parameters = dist.params, n.param.est = 2, method = method, 
        data.name = data.name, bad.obs = bad.obs)
    oldClass(ret.list) <- "estimate"
    ret.list
}
