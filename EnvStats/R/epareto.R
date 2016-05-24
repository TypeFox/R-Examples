epareto <-
function (x, method = "mle", plot.pos.con = 0.375) 
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
    if (n < 2 || any(x <= 0) || length(unique(x)) < 2) 
        stop(paste("'x' must contain at least 2 non-missing distinct values,", 
            "and all non-missing values of x must be positive."))
    if (!is.vector(plot.pos.con, mode = "numeric") || length(plot.pos.con) != 
        1 || plot.pos.con < 0 || plot.pos.con > 1) 
        stop("'plot.pos.con' must be a numeric scalar between 0 and 1")
    method <- match.arg(method, c("mle", "lse"))
    switch(method, lse = {
        F.hat <- ppoints(n, a = plot.pos.con)
        coef <- lm(log(1 - F.hat) ~ log(sort(x)))$coefficients
        shape <- -coef[2]
        location <- exp(coef[1]/shape)
        dist.params <- c(location, shape)
    }, mle = {
        location <- min(x)
        shape <- n/sum(log(x/location))
        dist.params <- c(location, shape)
    })
    names(dist.params) <- c("location", "shape")
    ret.list <- list(distribution = "Pareto", sample.size = n, 
        parameters = dist.params, method = method, data.name = data.name, 
        bad.obs = bad.obs)
    oldClass(ret.list) <- "estimate"
    ret.list
}
