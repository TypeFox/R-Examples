aveISFE = function (data, order, N, h, method, fmethod, lambda, mean, level,
    ...)
{
    n <- ncol(data$y)
    shortdata <- nextdata <- data
    m <- max(h)
    isfe <- matrix(NA, n - m - N + 1, m)
    for (i in N:(n - m)) {
        shortdata$y <- data$y[, 1:i]
        nextdata$y <- data$y[, i + (1:m)]
        fit <- ftsm(shortdata, order = order, method = method,
            lambda = lambda, mean = mean, level = level)
        fcast <- forecast(fit, h = m, ...)
        isfe[i - N + 1, ] <- MISE(nextdata, fcast$mean)$MISE
    }
    return(isfe[, h])
}
