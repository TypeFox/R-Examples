p0 <-
function (x, distr = NULL, mu = NULL, sigma = NULL) 
{
    if (is.null(distr)) {
        stop("Argument \"distr\" should be defined numeric with possible values 1,2, or 3")
    }
    else if (distr == 1) {
        a <- ifelse(is.null(mu), 0, mu)
        b <- ifelse(is.null(sigma), 1, sigma)
        p0 <- dnorm(x, mean = a, sd = b)
    }
    else if (distr == 2) {
        a <- ifelse(is.null(mu), 1, mu^2/sigma^2)
        b <- ifelse(is.null(sigma), 1, mu/sigma^2)
        p0 <- dgamma(x, shape = a, rate = b)
    }
    else if (distr == 3) {
        a <- ifelse(is.null(mu), 0.5, (1 - mu) * (mu/sigma)^2 - 
            mu)
        b <- ifelse(is.null(sigma), 1/sqrt(12), (mu * (1 - mu)/sigma^2 - 
            1) * (1 - mu))
        if (any(c(a, b) <= 0)) 
            stop(paste("\nNegative Beta parameters:\n a =", a, 
                ";\t b =", b))
        p0 <- dbeta(x, shape1 = a, shape2 = b)
    }
    else {
        stop("Argument \"distr\" should be defined numeric with possible values 1,2, or 3")
    }
    return(p0)
}
