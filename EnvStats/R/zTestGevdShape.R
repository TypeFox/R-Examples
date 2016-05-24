zTestGevdShape <-
function (x, pwme.method = "unbiased", plot.pos.cons = c(a = 0.35, 
    b = 0), alternative = "two.sided") 
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
    if (n < 3 || length(unique(x)) < 3) 
        stop(paste("'x' must contain at least 3 non-missing distinct values. ", 
            "This is not true for 'x' =", data.name))
    pwme.method <- match.arg(pwme.method, c("unbiased", "plotting.position"))
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    est.list <- egevd(x, method = "pwme", pwme.method = pwme.method, 
        plot.pos.cons = plot.pos.cons)
    shape.hat <- est.list$parameters["shape"]
    sd.shape.hat <- sqrt(0.5633/n)
    ret.list <- z.test.normal.approx(theta.hat = shape.hat, sd.theta.hat = sd.shape.hat, 
        hyp.theta = 0, alternative = alternative)
    names(ret.list$estimate) <- "shape"
    names(ret.list$null.value) <- "shape"
    ret.list$method <- "Z-test of shape=0 for GEVD"
    ret.list <- c(ret.list, list(estimation.method = est.list$method, 
        sample.size = n, data.name = data.name, bad.obs = bad.obs))
    oldClass(ret.list) <- "htest"
    ret.list
}
