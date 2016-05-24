`[.expectreg` <-
function (x, i) 
{
    stopifnot(i > 0 && i < length(effects(x)))
    lambda = NULL
    coefficients = NULL
    values = NULL
    covariates = NULL
    effects = NULL
    helper = NULL
    for (k in i) {
        lambda = c(lambda, x$lambda[[k]])
        coefficients = c(coefficients, x$coef[[k]])
        values = c(values, x$values[[k]])
        covariates = c(covariates, x$cov[[k]])
        effects = c(effects, x$effects[[k]])
        helper = c(helper, x$helper[[k]])
    }
    x$intercepts = x$intercepts[i]
    x$lambda = lambda
    x$coefficients = coefficients
    x$values = values
    x$covariates = covariates
    x$effects = effects
    x$helper = helper
    return(x)
}
