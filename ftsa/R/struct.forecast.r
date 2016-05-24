struct.forecast = function (x, h = 10, level = c(80, 95))
{
    fit1 <- StructTS(x, "level")
    fit2 <- StructTS(x, "trend")
    if (-2 * fit1$loglik < -2 * fit2$loglik + 2)
        fitStruct <- fit1
    else fitStruct <- fit2
    pred <- predict(fitStruct, n.ahead = h)
    nint <- length(level)
    lower <- matrix(NA, ncol = nint, nrow = length(pred$pred))
    upper <- lower
    for (i in 1:nint) {
        qq <- qnorm(0.5 * (1 + level[i]/100))
        lower[, i] <- pred$pred - qq * pred$se
        upper[, i] <- pred$pred + qq * pred$se
    }
    colnames(lower) = colnames(upper) = paste(level, "%", sep = "")
    fits <- rowSums(fitStruct$fitted)
    tsp(fits) <- tsp(x)
    return(structure(list(method = "Structural local linear",
        model = fitStruct, level = level, mean = pred$pred, var = pred$se^2,
        lower = lower, upper = upper, x = x, fitted = fits),
        class = "forecast"))
}
