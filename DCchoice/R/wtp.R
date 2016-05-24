wtp <- function(object, b = NULL, bid = NULL, dist = NULL)
{
    if (class(object) == "sbchoice" | class(object) == "dbchoice") {
        X <- object$covariates
        b <- object$coefficients
        bid <- object$bid
        dist <- object$dist
    } else {
        X <- object
    }

    coef <- b
    names(coef) <- NULL
    npar <- length(coef)
    b <- coef[npar]
    Xb <- sum(colMeans(X) * coef[-npar])

    if (dist == "log-logistic") {
        func <- function(x) plogis(-(Xb + b * log(x)), lower.tail = FALSE)
        medianWTP <- exp(-Xb/b)
        meanWTP <- ifelse(abs(b) > 1, integrate(func, 0, Inf, stop.on.error = FALSE)$value, Inf)
        trunc.meanWTP <- integrate(func, 0, exp(max(bid)), stop.on.error = FALSE)$value
        adj.trunc.meanWTP <- integrate(func, 0, exp(max(bid)), stop.on.error = FALSE)$value/plogis(-(Xb + b * max(bid)))
    } else if (dist == "log-normal") {
        func <- function(x) pnorm(-(Xb + b * log(x)), lower.tail = FALSE)
        medianWTP <- exp(-Xb/b)
        meanWTP <- integrate(func, 0, Inf, stop.on.error = FALSE)$value
        trunc.meanWTP <- integrate(func, 0, exp(max(bid)), stop.on.error = FALSE)$value
        adj.trunc.meanWTP <- integrate(func, 0, exp(max(bid)), stop.on.error = FALSE)$value/pnorm(-(Xb + b * max(bid)))
    } else if (dist == "logistic") {
        func <- function(x) plogis(-(Xb + b * x), lower.tail = FALSE)
        medianWTP <- -Xb/b
        meanWTP <- integrate(func, 0, Inf, stop.on.error = FALSE)$value
        trunc.meanWTP <- integrate(func, 0, max(bid), stop.on.error = FALSE)$value
        adj.trunc.meanWTP <- integrate(func, 0, max(bid), stop.on.error = FALSE)$value/plogis(-(Xb + b * max(bid)))
    } else if (dist == "normal") {
        func <- function(x) pnorm(-(Xb + b * x), lower.tail = FALSE)
        medianWTP <- -Xb/b
        meanWTP <- integrate(func, 0, Inf, stop.on.error = FALSE)$value
        trunc.meanWTP <- integrate(func, 0, max(bid), stop.on.error = FALSE)$value
        adj.trunc.meanWTP <- integrate(func, 0, max(bid), stop.on.error = FALSE)$value/pnorm(-(Xb + b * max(bid)))
    } else if (dist == "weibull") {
        func <- function(x) pweibull(exp(-Xb - b * log(x)), shape = 1, lower.tail = FALSE)
        medianWTP <- exp(-Xb/b) * (log(2))^(-1/b)
        meanWTP <- ifelse(abs(b) > 1, exp(-Xb/b) * gamma(1 - 1/b), Inf)
        trunc.meanWTP <- integrate(func, 0, exp(max(bid)), stop.on.error = FALSE)$value
        adj.trunc.meanWTP <- integrate(func, 0, exp(max(bid)),
                                       stop.on.error = FALSE)$value/pweibull(exp(-Xb - b * max(bid)), shape = 1)
    }

    output <- list(meanWTP = meanWTP, trunc.meanWTP = trunc.meanWTP, 
                   adj.trunc.meanWTP = adj.trunc.meanWTP, medianWTP = medianWTP)
    return(output)
}
