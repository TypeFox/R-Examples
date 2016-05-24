loglike.norm <-
function (x, mean, sd) 
{
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    if (!is.vector(mean, mode = "numeric") || length(mean) != 
        1 || !is.finite(mean)) 
        stop("'mean' must be a non-missing, finite numeric scalar")
    if (!is.vector(sd, mode = "numeric") || length(sd) != 1 || 
        !is.finite(sd) || sd < 0) 
        stop("'sd' must be a non-missing, finite, positive numeric scalar")
    data.name <- deparse(substitute(x))
    n <- length(x)
    if (any(is.na(x))) 
        statistic <- NA
    else {
        statistic <- (-n/2) * log(2 * pi) - n * log(sd) - (1/(2 * 
            sd^2)) * sum((x - mean)^2)
    }
    names(statistic) <- "Log-Likelihood"
    parameters <- c(mean, sd)
    names(parameters) <- c("mean", "sd")
    list(distribution = "Normal", dist.abb = "norm", distribution.parameters = parameters, 
        statistic = statistic, sample.size = n, data.name = data.name)
}
