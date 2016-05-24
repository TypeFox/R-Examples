pointwise <-
function (results.predict, coverage = 0.99, simultaneous = FALSE, 
    individual = FALSE) 
{
    fit <- results.predict$fit
    df <- results.predict$df
    if (is.null(df)) 
        stop(paste("The argument 'results.predict' must be", 
            "the result of calling the predict() function", "with se.fit=TRUE"))
    alpha <- 1 - coverage
    p <- results.predict$n.coefs
    if (simultaneous && is.null(p)) {
        stop(paste("You can only set simultaneous=TRUE", "when the argument results.predict has a", 
            "component called n.coefs (e.g., as for the results", 
            "of predict.lm)."))
    }
    if (!individual) {
        if (!simultaneous) {
            con <- qt(1 - alpha/2, df)
        }
        else {
            con <- sqrt(p * qf(coverage, p, df))
        }
        limits <- con * results.predict$se.fit
    }
    else {
        if (!simultaneous) {
            con <- qt(1 - alpha/2, df)
            limits <- con * sqrt(results.predict$residual.scale^2 + 
                results.predict$se.fit^2)
        }
        else {
            limits <- sqrt(p * qf(1 - alpha/2, p, df)) * results.predict$se.fit + 
                results.predict$residual.scale * qnorm(1 - alpha/2) * 
                  sqrt(df/qchisq(1 - alpha/2, df))
        }
    }
    list(upper = fit + limits, fit = fit, lower = fit - limits)
}
