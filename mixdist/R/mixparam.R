## last modified June 2002

mixparam <- function(mu, sigma, pi = NULL) 
{
    k <- length(mu)
    if (is.null(pi)) 
        pi <- rep(1/k, k)
    if (length(sigma) == 1) 
        sigma <- rep(sigma, k)
    else if (length(sigma) > k) {
        warning(paste("The length of sigma is greater than that of mu, the first ", 
            k, " sigmas was used", sep = ""))
        sigma <- sigma[1:k]
    }
    else if (length(sigma) < k) {
        warning("The length of sigma is less than that of mu, the first sigma was used")
        sigma <- rep(sigma[1], k)
    }
    if (k > 1) {
        if (sum(mu[-k] - mu[-1] > 0) > 0) 
            stop("Means must be in ascending order.")
        if (sum(mu[-k] - mu[-1] == 0 & sigma[-k] - sigma[-1] == 
            0) > 0) 
            stop("Sigmas must be in ascending order when means are equal.")
    }
    data.frame(pi = pi, mu = mu, sigma = sigma)
}
