DearBeggLoglik <- function(w, theta, sigma, y, u, hij, lam){
        
    n <- length(y)
    k <- 1 + floor(n / 2)
    eta <- sqrt(u ^ 2 + sigma ^ 2)
    lambda <- c(lam, rep(2, k - 2), 1)
    lambda[k] <- lambda[k] + (n %% 2)

    ## compute Ai
    Ai <- rep(NA, n)
    for (i in 1:n){Ai[i] <- sum(hij[i, ] * w)}

    ## compute log-likelihood
    LL <- sum(lambda * log(w)) + sum(dnorm(y, mean = theta, sd = eta, log = TRUE)) - sum(log(Ai))
        
    return(list("LL" = LL))    
}
