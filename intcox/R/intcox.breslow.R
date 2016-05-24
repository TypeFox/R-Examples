intcox.breslow <-
function (formula, data, covar)             # Breslow-Estimator
{
    lokal.cens <- data$cens/3
    formula.covar <- formula[[3]]
    formula <- as.formula(Surv(data$mix, lokal.cens) ~ .)
    formula[[3]] <- formula.covar
    fit <- coxph(formula, data)$coef        # Cox-regression coefficient
    e <- exp(t(as.matrix(fit)) %*% t(covar))
    event.num <- cumsum(data$cens == 3)
    hazard.rate0 <- NULL
    for (i in 1:max(event.num)) hazard.rate0 <- c(hazard.rate0,
        1/sum(e[i <= event.num]))
    cumhaz <- cumsum(hazard.rate0)
    breslow.ret <- list(cumhaz = cumhaz, fit = fit)
    return(breslow.ret)
}
