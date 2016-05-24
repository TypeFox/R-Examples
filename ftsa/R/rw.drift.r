rw.drift = function (x)
{
    x.diff <- diff(x)
    model <- summary(lm(x.diff ~ 1))
    drift <- model$coefficients[1, 1]
    fits <- ts(x - c(NA, model$residuals))
    tsp(fits) <- tsp(x)
    return(list(drift = drift, sec = model$coefficients[1, 2],
        see = model$sigma, fits = fits))
}
