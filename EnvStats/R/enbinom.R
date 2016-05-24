enbinom <-
function (x, size, method = "mle/mme") 
{
    if (!is.vector(x, mode = "numeric") || is.factor(x) || !is.vector(size, 
        mode = "numeric") || is.factor(size)) 
        stop("'x' and 'size' must be numeric vectors")
    data.name <- paste(deparse(substitute(x)), deparse(substitute(size)), 
        sep = ", ")
    lx <- length(x)
    lsize <- length(size)
    pass <- (lx == lsize) || (lsize == 1)
    if (!pass) 
        stop(paste("The length 'size' must be 1 or the same", 
            "length as 'x'."))
    if (lx > lsize) 
        size <- rep(size, lx)
    if ((bad.obs <- sum(!(all.ok <- is.finite(x) & is.finite(size)))) > 
        0) {
        is.not.finite.warning(x)
        is.not.finite.warning(size)
        x <- x[all.ok]
        size <- size[all.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' and/or 'size' removed."))
    }
    n <- length(x)
    if (n < 1) 
        stop("'x' and 'size' must contain at least one non-missing pair of values.")
    if (!all(x == trunc(x)) || any(x < 0) || !all(size == trunc(size)) || 
        any(size < 1)) 
        stop(paste("All values of 'x' must be non-negative integers,", 
            "and all values of 'size' must be positive integers."))
    x <- sum(x)
    size <- sum(size)
    method <- match.arg(method, c("mle/mme", "mvue"))
    dist.params <- switch(method, `mle/mme` = c(size = size, 
        prob = size/(size + x)), mvue = c(size = size, prob = (size - 
        1)/(size + x - 1)))
    ret.list <- list(distribution = "Negative Binomial", sample.size = n, 
        parameters = dist.params, n.param.est = 1, method = paste(method, 
            "for 'prob'"), data.name = data.name, bad.obs = bad.obs)
    oldClass(ret.list) <- "estimate"
    ret.list
}
