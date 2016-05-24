ppccGofTest <-
function (x, distribution = c("norm", "lnorm", "lnormAlt", "zmnorm", 
    "zmlnorm", "zmlnormAlt"), est.arg.list = NULL) 
{
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    data.name <- deparse(substitute(x))
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    distribution <- match.arg(distribution)
    if (any(distribution == c("lnorm", "lnormAlt")) && any(x <= 
        0)) 
        stop("All values of 'x' must be positive for a lognormal distribution")
    if (any(distribution == c("zmlnorm", "zmlnormAlt")) && any(x < 
        0)) 
        stop(paste("All values of 'x' must be non-negative for a", 
            "zero-modified lognormal distribution"))
    n <- switch(distribution, zmnorm = sum(x != 0), zmlnorm = , 
        zmlnormAlt = sum(x > 0), length(x))
    if (n < 5 || length(unique(x)) < 2) 
        stop(paste("'x' must contain at least 5 non-missing values,", 
            "and at least 2 distinct values. ", "This is not true for 'x' =", 
            data.name))
    if (n > 5000) 
        warning(paste("Too many observations.  This approximation only works", 
            "if the number of observations is between 5 and 5000"))
    gof.list <- do.call("sfGofTest", c(list(x = x, distribution = distribution), 
        est.arg.list))
    gof.list$method <- "PPCC GOF"
    gof.list$data.name <- data.name
    r <- sqrt(gof.list$statistic)
    names(r) <- "r"
    gof.list$statistic <- r
    gof.list
}
