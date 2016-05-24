ConvertWeibull <- function (model, conf.level = 0.95){

    level <- 1 - conf.level
    qa <- qnorm(1 - level / 2)
    Int.Only <- (nrow(summary(model)$table) == 2)
    sigma <- summary(model)$scale
    mu <- summary(model)$coef[1]
    k <- length(summary(model)$coef) - 1
    if (!Int.Only){gamma <- summary(model)$coef[2:(k + 1)]}

    lambda <- exp(-mu / sigma)
    alpha <- 1 / sigma
    tmp <- c(lambda, alpha)
    names(tmp) <- c("lambda", "gamma")
    if (!Int.Only){
        beta <- -gamma / sigma
        tmp <- c(lambda, alpha, beta)
        names(tmp) <- c("lambda", "gamma", names(summary(model)$coef[2:(k + 1)]))
    }

    var1 <- summary(model)$var
    var.mu <- diag(var1)[1]
    var.sigma <- var1[(k + 2), (k + 2)] * exp(2 * log(sigma))
    if (!Int.Only){
        var.gamma <- var1[2:(k + 1), 2:(k + 1)]
        if(k > 1){var.gamma <- diag(var.gamma)}
        se.gamma <- sqrt(var.gamma)
    }

    cov.mu.sigma <- var1[(k + 2), 1] * sigma
    var.alpha <- var.sigma / (sigma ^ 4)
    var.lambda <- exp(-2 * mu / sigma) * ((var.mu / (sigma ^ 2)) -
        ((2 * mu / (sigma ^ 3)) * cov.mu.sigma) + (((mu ^ 2)/(sigma ^ 4)) * var.sigma))

    var <- c(sqrt(var.lambda), sqrt(var.alpha))

    if (!Int.Only){
        cov.gamma.sigma <- var1[2:(k + 1), (k + 2)] * sigma
        var.beta <- (1 / (sigma ^ 2)) * (var.gamma - (2 * gamma / sigma) *
            (cov.gamma.sigma) + (((gamma/sigma) ^ 2) * var.sigma))
        se.beta <- sqrt(var.beta)
        var <- c(sqrt(var.lambda), sqrt(var.alpha), se.beta)
        HR <- cbind(HR = exp(beta), LB = exp(beta - qa * se.beta), UB = exp(beta + qa * se.beta))
        rownames(HR) <- names(summary(model)$coef[2:(k + 1)])
        ETR <- cbind(ETR = exp(gamma), LB = exp(gamma - qa * se.gamma), UB = exp(gamma + qa * se.gamma))
        rownames(HR) <- names(summary(model)$coef[2:(k + 1)])
    }

    tmp1 <- rbind(tmp, var)
    rownames(tmp1) <- c("Estimate", "SE")
    ret <- list(vars = t(tmp1))
    if (!Int.Only){
        ret$HR <- HR
        ret$ETR <- ETR
    }

    return(ret)
}


