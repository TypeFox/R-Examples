egeom <-
function (x, method = "mle/mme") 
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
    if (n < 1) 
        stop("'x' must contain at least one non-missing value.")
    if (!all(x == trunc(x)) || any(x < 0)) 
        stop("All values of 'x' must be non-negative integers.")
    x <- sum(x)
    method <- match.arg(method, c("mle/mme", "mvue"))
    if (method == "mvue" && n < 2) 
        stop(paste("You must have at least two non-missing", 
            "values in 'x' to compute the 'mvue' of 'prob'."))
    dist.params <- switch(method, `mle/mme` = c(prob = n/(n + 
        x)), mvue = c(prob = (n - 1)/(n + x - 1)))
    ret.list <- list(distribution = "Geometric", sample.size = n, 
        parameters = dist.params, n.param.est = 1, method = method, 
        data.name = data.name, bad.obs = bad.obs)
    oldClass(ret.list) <- "estimate"
    ret.list
}
