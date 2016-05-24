pegelsna = function (x, h = 10, level = c(80, 95), upper = c(0.3, 0.2, 0.99),
    model = "AZN")
{
    MSE1 <- function(alpha, x) {
        if (alpha > upper[1])
            return(1e+09)
        if (alpha < 0.01)
            return(1e+09)
        fit <- Arima(x, order = c(0, 1, 1), fixed = alpha - 1)
        return(fit$sigma2)
    }
    MSE2 <- function(alpha, x) {
        if (alpha[1] > upper[1] | alpha[2] > upper[2])
            return(1e+09)
        if (min(alpha) < 0.01)
            return(1e+09)
        theta1 <- alpha[1] + alpha[1] * alpha[2] - 2
        theta2 <- 1 - alpha[1]
        fit <- Arima(x, order = c(0, 2, 2), fixed = c(theta1,
                     theta2))
        return(fit$sigma2)
    }
    MSE3 <- function(alpha, x) {
        if (alpha[1] > upper[1] | alpha[2] > upper[2])
            return(1e+09)
        if (min(alpha) < 0.01)
            return(1e+09)
        if (alpha[3] > upper[3])
            return(1e+09)
        if (alpha[3] < alpha[2])
            return(1e+09)
        theta1 <- alpha[1] + alpha[1] * alpha[2] - 1 - alpha[3]
        theta2 <- (1 - alpha[1]) * alpha[3]
        phi1 <- alpha[3]
        fit <- Arima(x, order = c(1, 1, 2), fixed = c(phi1, theta1,
                     theta2))
        return(fit$sigma2)
    }
    n <- length(x)
    if (model == "ANN") {
        fit1 <- nlm(MSE1, upper[1]/2, x = x)
        fitarima <- Arima(x, order = c(0, 1, 1), fixed = fit1$estimate -
            1)
        model = structure(list(alpha = fit1$estimate, beta = 0,
            gamma = 0, phi = 1, initstate = NA, rmse = sqrt(fitarima$sigma2)),
            class = "ets")
        method <- "Robust SES"
    }
    else if (model == "AAN") {
        fit2 <- nlm(MSE2, upper[1:2]/2, x = x)
        theta1 <- fit2$estimate[1] + fit2$estimate[1] * fit2$estimate[2] -
            2
        theta2 <- 1 - fit2$estimate[1]
        fitarima <- Arima(x, order = c(0, 2, 2), fixed = c(theta1,
            theta2))
        model = structure(list(alpha = fit2$estimate[1], beta = fit2$estimate[2],
            gamma = 0, phi = 1, initstate = NA, rmse = sqrt(fitarima$sigma2)),
            class = "ets")
        method <- "Robust Holt's"
    }
    else if (model == "ADN") {
        fit3 <- nlm(MSE3, upper * c(0.5, 0.5, 0.95), x = x)
        theta1 <- fit3$estimate[1] + fit3$estimate[1] * fit3$estimate[2] -
                  1 - fit3$estimate[3]
        theta2 <- (1 - fit3$estimate[1]) * fit3$estimate[3]
        phi1 <- fit3$estimate[3]
        fitarima <- Arima(x, order = c(1, 1, 2), fixed = c(phi1,
                          theta1, theta2))
        model = structure(list(alpha = fit3$estimate[1], beta = fit3$estimate[2],
                          gamma = 0, phi = fit3$estimate[3], initstate = NA,
                          rmse = sqrt(fitarima$sigma2)), class = "ets")
        method <- "Robust Damped Holt's"
    }
    else if (model == "AZN") {
        fit1 <- nlm(MSE1, upper[1]/2, x = x)
        fit2 <- nlm(MSE2, upper[1:2]/2, x = x)
        fit3 <- nlm(MSE3, upper * c(0.5, 0.5, 0.95), x = x)
        n <- length(x)
        aic <- c(n * log(fit1$minimum) + 2, n * log(fit2$minimum) +
                 4, n * log(fit3$minimum) + 6)
        best <- order(aic)[1]
        if (best == 1) {
            fitarima <- Arima(x, order = c(0, 1, 1), fixed = fit1$estimate -
                              1)
            model = structure(list(alpha = fit1$estimate, beta = 0,
                gamma = 0, phi = 1, initstate = NA, rmse = sqrt(fitarima$sigma2)),
                class = "ets")
            method <- "Robust SES"
        }
        else if (best == 2) {
            theta1 <- fit2$estimate[1] + fit2$estimate[1] * fit2$estimate[2] -
                2
            theta2 <- 1 - fit2$estimate[1]
            fitarima <- Arima(x, order = c(0, 2, 2), fixed = c(theta1,
                theta2))
            model = structure(list(alpha = fit2$estimate[1],
                beta = fit2$estimate[2], gamma = 0, phi = 1,
                initstate = NA, rmse = sqrt(fitarima$sigma2)),
                class = "ets")
            method <- "Robust Holt's"
        }
        else if (best == 3) {
            theta1 <- fit3$estimate[1] + fit3$estimate[1] * fit3$estimate[2] -
                1 - fit3$estimate[3]
            theta2 <- (1 - fit3$estimate[1]) * fit3$estimate[3]
            phi1 <- fit3$estimate[3]
            fitarima <- Arima(x, order = c(1, 1, 2), fixed = c(phi1,
                theta1, theta2))
            model = structure(list(alpha = fit3$estimate[1],
                beta = fit3$estimate[2], gamma = 0, phi = fit3$estimate[3],
                initstate = NA, rmse = sqrt(fitarima$sigma2)),
                class = "ets")
            method <- "Robust Damped Holt's"
        }
    }
    else stop("Unknown model")
    pred <- predict(fitarima, n.ahead = h)
    nint <- length(level)
    lower <- matrix(NA, ncol = nint, nrow = length(pred$pred))
    upper <- lower
    for (i in 1:nint) {
        qq <- qnorm(0.5 * (1 + level[i]/100))
        lower[, i] <- pred$pred - qq * pred$se
        upper[, i] <- pred$pred + qq * pred$se
    }
    colnames(lower) = colnames(upper) = paste(level, "%", sep = "")
    f = frequency(x)
    return(structure(list(method = method, model = model, level = level,
        mean = pred$pred, var = pred$se^2, lower = lower, upper = upper,
        x = x, fitted = fitted(fitarima)), class = "forecast"))
}