dk <-
function (x, distr = NULL, mu = NULL, sigma = NULL) 
{
    if (is.null(distr)) {
        stop("Argument \"distr\" should be defined numeric with possible values 1,2,3,4 or 5")
    }
    else if (distr == 1) {
        a <- ifelse(is.null(mu), 0, mu)
        b <- ifelse(is.null(sigma), 1, sigma)
        dk <- dnorm(x, mean = a, sd = b)
    }
    else if (distr == 4) {
        a <- ifelse(is.null(mu), 0, mu)
        b <- ifelse(is.null(sigma), 1/sqrt(2), sigma/sqrt(2))
        dk <- exp(-abs(x - a)/b)/(2 * b)
    }
    else if (distr == 5) {
        a <- ifelse(is.null(mu), exp(1/2), log(mu/sqrt(1 + (sigma/mu)^2)))
        b <- ifelse(is.null(sigma), exp(1) * (exp(1) - 1), sqrt(log(1 + 
            (sigma/mu)^2)))
        dk <- dlnorm(x, meanlog = a, sdlog = b)
    }
    else if (distr == 2) {
        a <- ifelse(is.null(mu), 1, mu^2/sigma^2)
        b <- ifelse(is.null(sigma), 1, mu/sigma^2)
        dk <- dgamma(x, shape = a, rate = b)
    }
    else if (distr == 3) {
        a <- ifelse(is.null(mu), 0.5, (1 - mu) * (mu/sigma)^2 - 
            mu)
        b <- ifelse(is.null(sigma), 1/sqrt(12), (mu * (1 - mu)/sigma^2 - 
            1) * (1 - mu))
        if (any(c(a, b) <= 0)) 
            stop(paste("\nNegative Beta parameters:\n a =", a, 
                ";\t b =", b))
        dk <- dbeta(x, shape1 = a, shape2 = b)
    }
    else {
        stop("Argument \"distr\" should be defined numeric with possible values 1,2,3,4 or 5")
    }
    return(dk)
}
