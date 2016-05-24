sfGeneralGofTest <-
function (x, distribution, est.arg.list) 
{
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    data.name <- deparse(substitute(x))
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    est.fcn <- paste("e", distribution, sep = "")
    est.list <- do.call(est.fcn, c(list(x = x), est.arg.list))
    params <- est.list$parameters
    Z <- do.call(paste("p", distribution, sep = ""), c(list(q = x), 
        as.list(params)))
    Y <- qnorm(Z)
    ret.list <- sfGofTest(Y)
    ret.list$data <- x
    ret.list$data.name <- data.name
    ret.list$bad.obs <- bad.obs
    ret.list$dist.abb <- distribution
    ret.list$distribution <- EnvStats::Distribution.df[distribution, 
        "Name"]
    ret.list$distribution.parameters <- params
    ret.list$n.param.est <- length(params)
    ret.list$estimation.method <- est.list$method
    ret.list$alternative <- paste("True cdf does not equal the\n", 
        space(33), ret.list$distribution, " Distribution.", sep = "")
    ret.list$method <- "Shapiro-Francia GOF Based on Chen & Balakrishnan (1995)"
    ret.list
}
