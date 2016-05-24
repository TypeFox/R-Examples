ebeta <-
function (x, method = "mle") 
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
    if (n < 2 || min(x) < 0 || max(x) > 1 || length(unique(x)) < 
        2) 
        stop(paste("'x' must contain at least 2 non-missing distinct values,", 
            "and all non-missing values of 'x' must be between 0 and 1."))
    m <- mean(x)
    term <- ((m * (1 - m))/(((n - 1)/n) * var(x))) - 1
    shape1 <- m * term
    shape2 <- (1 - m) * term
    method <- match.arg(method, c("mle", "mme", "mmue"))
    switch(method, mme = {
        dist.params <- c(shape1 = shape1, shape2 = shape2)
    }, mmue = {
        term <- ((m * (1 - m))/var(x)) - 1
        shape1 <- m * term
        shape2 <- (1 - m) * term
        dist.params <- c(shape1 = shape1, shape2 = shape2)
    }, mle = {
        fcn <- function(theta, mlx, ml1mx) {
            term1 <- digamma(theta[1]) - digamma(sum(theta)) - 
                mlx
            term2 <- digamma(theta[2]) - digamma(sum(theta)) - 
                ml1mx
            (term1^2) + (term2^2)
        }
        dist.params <- nlminb(start = c(shape1, shape2), objective = fcn, 
            mlx = mean(log(x)), ml1mx = mean(log(1 - x)), lower = .Machine$double.eps)$par
        names(dist.params) <- c("shape1", "shape2")
    })
    ret.list <- list(distribution = "Beta", sample.size = n, 
        parameters = dist.params, n.param.est = 2, method = method, 
        data.name = data.name, bad.obs = bad.obs)
    oldClass(ret.list) <- "estimate"
    ret.list
}
