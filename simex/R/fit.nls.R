fit.nls <-
function (lambda, p.names, estimates)
{
    extrapolation <- list()
    # lambda values for the interpolation
    lambdastar <- c(0, max(lambda) / 2, max(lambda))
    for (d in p.names) {
        # calculating starting values
        # quadratic interpolation
        extrapolation.quad <- lm(estimates[, d] ~ lambda + I(lambda^2))
        # interpolation for the values of lambdastar
        a.nls <- predict(extrapolation.quad,
            newdata = data.frame(lambda = lambdastar))
        # analytic 3-point fit => good starting values
        gamma.est.3 <- ((a.nls[2] - a.nls[3]) * lambdastar[3] *
            (lambdastar[2] - lambdastar[1]) - lambdastar[1] *
            (a.nls[1] - a.nls[2]) * (lambdastar[3] - lambdastar[2])) /
            ((a.nls[1] - a.nls[2]) * (lambdastar[3] - lambdastar[2]) -
             (a.nls[2] - a.nls[3]) * (lambdastar[2] - lambdastar[1]))
        gamma.est.2 <- ((a.nls[2] - a.nls[3]) * (gamma.est.3 + lambdastar[2]) *
            (gamma.est.3 + lambdastar[3])) / (lambdastar[3] - lambdastar[2])
        gamma.est.1 <- a.nls[1] - (gamma.est.2 / (gamma.est.3 + lambdastar[1]))
        # fitting the nls-model for the various coefficients
        extrapolation[[d]] <-
            nls(estimates[, d] ~ gamma.1 + gamma.2 / (gamma.3 + lambda),
                start = list(gamma.1 = gamma.est.1, gamma.2 = gamma.est.2,
                gamma.3 = gamma.est.3))
    }
    return(extrapolation)
}

