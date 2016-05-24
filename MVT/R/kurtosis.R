kurtosis <- function(x, center, cov)
{   # coefficient of multivariate kurtosis (Mardia, 1970)
    distances <- mahalanobis(x, center, cov)
    b <- mean(distances^2)
    # estimation of the kurtosis parameter for the t-distribution
    p <- nrow(cov)
    k <- b / (p * (p + 2)) - 1
    if (k < 0)
        k <- NaN
    eta <- .5 * k / (1 + 2 * k)
    list(kurtosis = b, kappa = k, eta = eta)
}
