CIsim <-
function (samples = 100, n = 30, parameter = 0.5, sigma = 1,
    conf.level = 0.95, type = "Mean")
{
    Adkblue <- "#0080FF"
    Aorange <- "#FF4C0C"
    alpha <- 1 - conf.level
    CL <- conf.level * 100
    N <- samples
    choices <- c("Mean", "Var", "Pi")
    alt <- pmatch(type, choices)
    type <- choices[alt]
    if (length(type) > 1 || is.na(type))
        stop("alternative must be one \"Mean\", \"Var\", \"Pi\"")
    if (type == "Pi" && (parameter <= 0 | parameter >= 1))
        stop("Value for Pi (parameter) must be between 0 and 1.")
    if (N <= 0 || n <= 0)
        stop("Number of random CIs (samples) and sample size (n) must both be at least 1")
    if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
        conf.level <= 0 || conf.level >= 1))
        stop("'conf.level' must be a single number between 0 and 1")
    if (sigma <= 0 && (type == "Var" || type == "Mean"))
        stop("Variance must be a positive value")
    if (type == "Mean") {
        junk <- rnorm(N * n, parameter, sigma)
        jmat <- matrix(junk, N, n)
        xbar <- apply(jmat, 1, mean)
        ll <- xbar - qnorm(1 - alpha/2) * sigma/sqrt(n)
        ul <- xbar + qnorm(1 - alpha/2) * sigma/sqrt(n)
        notin <- sum((ll > parameter) + (ul < parameter))
        percentage <- round((notin/N) * 100, 2)
        plot(ll, type = "n", ylim = c(min(ll), max(ul)), xlab = " ",
            ylab = " ")
        title(sub = bquote(paste("Note: ", .(percentage), "% of the random confidence intervals do not contain ",
            mu, "=", .(parameter))))
        title(main = bquote(paste(.(N), " random ", .(CL), "% confidence intervals where ",
            mu, " = ", .(parameter))))
        for (i in 1:N) {
            low <- ll[i]
            high <- ul[i]
            if (low < parameter & high > parameter) {
                segments(i, low, i, high)
            }
            else if (low > parameter & high > parameter) {
                segments(i, low, i, high, col = Aorange, lwd = 5)
            }
            else {
                segments(i, low, i, high, col = Adkblue, lwd = 5)
            }
        }
        abline(h = parameter)
        cat(percentage, "% of the random confidence intervals do not contain Mu =",
            parameter, ".", "\n")
    }
    else if (type == "Var") {
        junk <- rnorm(N * n, parameter, sigma)
        jmat <- matrix(junk, N, n)
        s2 <- apply(jmat, 1, var)
        ll <- ((n - 1) * s2)/qchisq(1 - alpha/2, (n - 1))
        ul <- ((n - 1) * s2)/qchisq(alpha/2, (n - 1))
        variance <- sigma^2
        notin <- sum((ll > variance) + (ul < variance))
        percentage <- round((notin/samples) * 100, 2)
        plot(ll, type = "n", ylim = c(min(ll), max(ul)), xlab = " ",
            ylab = " ")
        title(sub = bquote(paste("Note: ", .(percentage), "% of the random confidence intervals do not contain ",
            sigma^2, "=", .(variance))))
        title(main = bquote(paste(.(N), " random ", .(CL), "% confidence intervals where ",
            sigma^2, " = ", .(variance))))
        for (i in 1:N) {
            low <- ll[i]
            high <- ul[i]
            if (low < variance & high > variance) {
                segments(i, low, i, high)
            }
            else if (low > variance & high > variance) {
                segments(i, low, i, high, col = Aorange, lwd = 5)
            }
            else {
                segments(i, low, i, high, col = Adkblue, lwd = 5)
            }
        }
        abline(h = variance)
        cat(percentage, "% of the random confidence intervals do not contain Var =",
            sigma^2, ".", "\n")
    }
    else if (type == "Pi") {
        X <- rbinom(samples, n, parameter)
        p <- X/n
        ll <- p - qnorm(1 - alpha/2) * sqrt((p * (1 - p))/n)
        ul <- p + qnorm(1 - alpha/2) * sqrt((p * (1 - p))/n)
        notin <- sum((ll > parameter) + (ul < parameter))
        percentage <- round((notin/samples) * 100, 2)
        plot(ll, type = "n", ylim = c(min(ll), max(ul)), xlab = " ",
            ylab = " ")
        title(sub = bquote(paste("Note: ", .(percentage), "% of the random confidence intervals do not contain ",
            pi, "=", .(parameter))))
        title(main = bquote(paste(.(N), " random ", .(CL), "% confidence intervals where ",
            pi, "=", .(parameter))))
        for (i in 1:N) {
            low <- ll[i]
            high <- ul[i]
            if (low < parameter & high > parameter) {
                segments(i, low, i, high)
            }
            else if (low > parameter & high > parameter) {
                segments(i, low, i, high, col = Aorange, lwd = 5)
            }
            else {
                segments(i, low, i, high, col = Adkblue, lwd = 5)
            }
        }
        abline(h = parameter)
        cat(percentage, "% of the random confidence intervals do not contain Pi =",
            parameter, ".", "\n")
    }
}

